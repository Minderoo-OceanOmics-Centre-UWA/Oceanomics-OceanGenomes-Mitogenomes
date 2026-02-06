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
                sample_list.collect { sample_record ->
                    // sample_record[0] is a meta map like [id:OG1341] so we want to destructure the map.
                    def raw_meta  = sample_record[0]
                    def sample_id = (raw_meta instanceof Map) ? raw_meta.id : raw_meta
                    
                    def meta = (raw_meta instanceof Map) ? raw_meta : [ id: raw_meta ]
                    meta = meta + [ sequencing_type: sample_record[3] ]
                    def fastq_1 = sample_record[1]
                    def fastq_2 = sample_record[2]
                    def single_end = (!fastq_2 || fastq_2.toString() == '[]')
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

                    if (single_end) {
                        return [ meta.id, meta, [ fastq_1 ] ]
                    }

                    return [ meta.id, meta, [ fastq_1, fastq_2 ] ]
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
            .view()

    // DEBUG
    ch_samplesheet.view { meta, _reads -> "Sample: ${meta.id}, Type: ${meta.sequencing_type}, Meta: ${meta}" }
    
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
