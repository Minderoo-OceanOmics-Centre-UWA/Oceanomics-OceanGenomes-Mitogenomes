//
// Subworkflow with functionality specific to the nf-core/oceangenomesmitogenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_SAMPLESHEET_ENRICHED } from '../../../modules/local/create_samplesheet_enriched'
include { samplesheetToList         } from 'plugin/nf-schema'

// Resolve the NCBI mitochondrial genetic code (translation table) for a sample
// from its taxonomic class. Mitochondrial codes differ by lineage, and the wrong
// code mistranslates CDS in MITOS2 annotation and QC protein translation:
//   * Cnidaria (corals, anemones, jellyfish, hydroids) use the Coelenterate code (4),
//   * echinoderms and flatworms use the Echinoderm/Flatworm code (9),
//   * any other invertebrate falls back to the Coelenterate code (4): in practice
//     the only invertebrates routed through this pipeline are cnidarians, so an
//     unrecognised invertebrate class is far more likely to be one of those than
//     a code-5 invertebrate. (Add an explicit class above if that ever changes.)
//   * vertebrates (and anything unresolved) fall back to defaultCode (vertebrate, 2).
def mitoGeneticCode(taxClass, isInvert, defaultCode) {
    def c = (taxClass ?: '').toString().trim().toLowerCase()
    if (c in ['anthozoa', 'hydrozoa', 'scyphozoa', 'cubozoa', 'staurozoa', 'myxozoa', 'polypodiozoa']) {
        return 4
    }
    if (c in ['asteroidea', 'ophiuroidea', 'echinoidea', 'holothuroidea', 'crinoidea',
              'rhabditophora', 'trematoda', 'cestoda', 'monogenea', 'turbellaria']) {
        return 9
    }
    if (isInvert) {
        return 4
    }
    return defaultCode
}

// Identify a read file that came from the "unassigned" bin of HiFi barcode
// demultiplexing (e.g. *.hifi_reads.unassigned.filt.fastq.gz). These reads failed
// demux and can belong to any specimen on the SMRT cell, so they must never be
// assembled into a sample even if a row for one lands on the samplesheet.
def isUnassignedReadFile(readPath) {
    if (!readPath) return false
    def base = readPath.toString().tokenize('/').last().toLowerCase()
    return base.contains('unassigned')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_SAMPLESHEET {

    take:
    input               //  string: Path to input samplesheet
    input_dir           //  string: Path to input samplesheet
    samplesheet_name 
    
    main:

   // Handle different input types
    if (input && input_dir) {
        error "Please specify either --input (samplesheet) OR --input_dir (directory), not both"
    }
        
    else if (input_dir) {
        if (!params.sql_config) {
            error "No --sql_config provided. When using --input_dir, please provide --input (samplesheet) instead."
        }
        // Resolve the glob pattern to actual files
        input_files_ch = Channel.fromPath(input_dir, checkIfExists: true)
            .collect()


        // Create enriched samplesheet from directory
        CREATE_SAMPLESHEET_ENRICHED(
            input_files_ch,
            "samplesheet.csv",
            params.sql_config
        )

        samplesheet_ch = CREATE_SAMPLESHEET_ENRICHED.out.samplesheet
    }
    else if (input) {
        samplesheet_ch = Channel.fromPath(input, checkIfExists: true)
    }
    else {
        error "No input provided: specify --input or --input_dir or provide reads downstream"
    }

    // Parse samplesheet into channel of [meta, fastq_files]
    def ch_samplesheet = samplesheet_ch
            .map { samplesheet_file ->
                samplesheetToList(samplesheet_file, "${projectDir}/assets/schema_input.json")
            }
            .flatMap { sample_list ->
                // Convert each sample record to the expected format
                // (findResults drops any record returned as null, e.g. unassigned reads)
                sample_list.findResults { sample_record ->
                    // sample_record[0] is a meta map like [id:OG1341] so we want to destructure the map.
                    def raw_meta  = sample_record[0]
                    def sample_id = (raw_meta instanceof Map) ? raw_meta.id : raw_meta
                    
                    def meta = (raw_meta instanceof Map) ? raw_meta : [ id: raw_meta ]

                    // nf-schema fills a column the samplesheet does not carry with an
                    // EMPTY LIST rather than omitting the key, so adding an optional
                    // column to schema_input.json silently changes the shape of every
                    // meta map -- even for sheets that never gained the column.
                    //
                    // That is not cosmetic. Nextflow's task hash ignores an empty
                    // collection, so a resumed run still hits the old cache and hands
                    // back the OLD meta, which then fails `equals` against the live one.
                    // Every join(..., by: 0) keyed on the whole meta map then matches
                    // nothing and drops the sample without an error or a warning. Adding
                    // `family` / `order` did exactly that: it emptied the MitoHiFi
                    // reference join for all 70 samples whose reads had not also come
                    // back from cache.
                    //
                    // Keep meta shape a function of what the samplesheet actually
                    // carries. Only the empty placeholders go; false / 0 / '' are real
                    // values and stay. Removing them changes no task hash, so cached
                    // work stays valid.
                    meta = meta.findAll { _key, value -> !(value instanceof Collection) || !value.isEmpty() }

                    meta = meta + [ sequencing_type: sample_record[3] ]
                    def fastq_1 = sample_record[1]
                    def fastq_2 = sample_record[2]
                    def single_end = (!fastq_2 || fastq_2.toString() == '[]')

                    // Defensive guard: never assemble from "unassigned" read files.
                    // These are reads that failed HiFi barcode demultiplexing and may
                    // belong to other specimens on the same SMRT cell. Drop such rows
                    // even if they slip onto the samplesheet (checked on both mates so
                    // a paired unassigned row is caught too). Toggle off with
                    // --exclude_unassigned_reads false.
                    if (params.exclude_unassigned_reads != false &&
                        (isUnassignedReadFile(fastq_1) || isUnassignedReadFile(fastq_2))) {
                        log.warn "Excluding unassigned read file(s) for ${sample_id} (${meta.sequencing_type}): ${fastq_1}${fastq_2 && fastq_2.toString() != '[]' ? ", ${fastq_2}" : ''}"
                        return null
                    }
                    def meta_single_end = meta.single_end
                    if (meta.date != null) {
                        meta = meta + [ date: meta.date.toString() ]
                    }

                    if (meta_single_end instanceof String) {
                        meta_single_end = meta_single_end.toLowerCase() == 'true'
                    }

                    if (meta_single_end == null || meta_single_end.toString().trim() == '') {
                        meta = meta + [ single_end: single_end ]
                    } else {
                        meta = meta + [ single_end: meta_single_end ]
                    }

                    def meta_invertebrates = meta.invertebrates

                    if (meta_invertebrates instanceof String) {
                        meta_invertebrates = meta_invertebrates.toLowerCase() == 'true'
                    }

                    if (meta_invertebrates == null || meta_invertebrates.toString().trim() == '') {
                        meta = meta + [ invertebrates: false ]
                    } else {
                        meta = meta + [ invertebrates: meta_invertebrates ]
                    }

                    // Derive the per-sample mitochondrial genetic code from the
                    // resolved taxonomic class so MITOS2 / QC translation use the
                    // correct code (e.g. Cnidaria -> 4) instead of a one-size
                    // global table. --translation_table sets the vertebrate/default.
                    def mt_genetic_code = mitoGeneticCode(meta.class, meta.invertebrates, (params.translation_table ?: 2) as int)
                    meta = meta + [ genetic_code: mt_genetic_code ]

                    // Group by sample id + sequencing type + date so that single-end
                    // (e.g. hifi bc + unassigned) and paired-end (e.g. hic lanes) reads
                    // for the same specimen don't collide on a shared cleaned id.
                    // Reads sharing this key are merged into one assembly downstream.
                    def group_key = [ meta.id, meta.sequencing_type, meta.date ]

                    if (single_end) {
                        return [ group_key, meta, [ fastq_1 ] ]
                    }

                    return [ group_key, meta, [ fastq_1, fastq_2 ] ]
                }
            }
            .groupTuple()
            .map { samplesheet ->
                return validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            // .view()

    // DEBUG
    // ch_samplesheet.view { meta, _reads -> "Sample: ${meta.id}, Type: ${meta.sequencing_type}, Meta: ${meta}" }
    
    // Branch by sequencing type
    ch_samplesheet.branch { meta, reads ->
        hifi: meta.sequencing_type == 'hifi'
            return [meta, reads]
        ilmn: meta.sequencing_type == 'ilmn'
            return [meta, reads]
        hic: meta.sequencing_type == 'hic'
            return [meta, reads]
    }.set { branched_samples }

    // Channels already contain enriched meta from the samplesheet
    ch_hifi_with_files = branched_samples.hifi
    ch_ilmn_with_meta  = branched_samples.ilmn
    ch_hic_with_files  = branched_samples.hic

    // Combine all processed data
    ch_getorg = ch_ilmn_with_meta.mix(ch_hic_with_files)

    emit:
    samplesheet = samplesheet_ch
    getorg_input = ch_getorg
    mitohifi_input = ch_hifi_with_files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (sample_id, metas, fastqs) = input

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
