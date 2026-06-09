/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../../nf-core/utils_nfcore_pipeline'

// Mitogenome assembly
include { GETORGANELLE_CONFIG       } from '../../../../modules/nf-core/getorganelle/config'
include { CAT_FASTQ                 } from '../../../../modules/nf-core/cat/fastq'
include { GETORGANELLE_FROMREADS    } from '../../../../modules/nf-core/getorganelle/fromreads'
include { GETORGANELLE_RESEED       } from '../../../../modules/local/getorganelle/reseed'
include { MITOHIFI_FINDMITOREFERENCE } from '../../../../modules/nf-core/mitohifi/findmitoreference'
include { PUSH_MTDNA_ASSM_RESULTS   } from '../../../../modules/local/upload_results/mtdna'

// Decide whether a first-pass GetOrganelle result warrants a reseed attempt.
// Triggered when the assembly is effectively absent (empty placeholder FASTA)
// or the log shows the classic "seed too divergent" markers and the run did NOT
// reach a circular genome. See get_org.log.txt markers observed in practice.
def needsReseed(fasta, logFile) {
    // No assembled contig at all -> always worth a reseed.
    try { if (fasta == null || fasta.size() == 0) return true } catch (ignored) { return true }
    def txt
    try { txt = logFile.text } catch (ignored) { return false }
    // A circular genome is a success regardless of mid-run warnings.
    if (txt =~ /Result status of .*: circular genome/) return false
    def markers = [
        'No sequence hit our LabelDatabase!',
        'unreasonable seed/parameter choices',
        'no target organelle contigs found',
        'insufficient number of rounds',
    ]
    return markers.any { txt.contains(it) }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ASSEMBLY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ASSEMBLY_GETORG {

    take:
    fastp_reads // tuple val(meta), path(fastp)
    organelle_type // "animal_mt" passed to download the correct congig files
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_summary_files = Channel.empty()

     
    //
    // MODULE: Set up for running GetOrganelle
    //

    GETORGANELLE_CONFIG (
        organelle_type // val(organelle_type)
    )

    // Split fastp_reads based on whether concatenation is needed
    fastp_reads_split = fastp_reads.branch { meta, reads ->
        def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
        def needsConcatenation = meta.single_end ? readList.size > 1 : readList.size > 2
        
        needs_concat: needsConcatenation
            return [meta, reads]
        no_concat: !needsConcatenation
            return [meta, reads]
    }
    //
    // MODULE: Concatenate fastq reads where there are multiple fastq R1 and R2 files
    //

    CAT_FASTQ (
        fastp_reads_split.needs_concat
    )
    
    // Combine the results
    final_reads = fastp_reads_split.no_concat.mix(CAT_FASTQ.out.reads)

    //
    // Get the GetOrganelle version to be used in the naming. Ensures the proper version is always used.
    //

    version_ch = GETORGANELLE_CONFIG.out.versions
    .map { versions_file ->
        def content = versions_file.text
        def pattern = /getorganelle:\s*([^\s\n\r]+)/
        def matcher = content =~ pattern
        return matcher ? matcher[0][1].replaceAll(/["']/, '') : "unknown"
    }

    //
    // Embed the assembly prefix into meta.
    //

    fastp_with_mt_assembly_prefix = final_reads.combine(version_ch)
    .map { meta, fasta, version ->  // Destructure all 3 elements correctly
        def version_stripped = version.replaceAll('\\.', '')
        def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.getorg${version_stripped}"
        def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
        [meta_ext, fasta]  // Return only what you want
    }


    //
    // Combine fastp files with the GetOrganelle database
    //   .first() makes the single-item db channel a reusable value channel so it
    //   can be combined into both the first pass and the reseed pass.
    //

    ch_db = GETORGANELLE_CONFIG.out.db.first()
    combined_input = fastp_with_mt_assembly_prefix.combine(ch_db)

    //
    // MODULE: Run assembly using GetOrganelle from reads (first pass, generic seed)
    //

    GETORGANELLE_FROMREADS (
        combined_input // tuple val(meta), path(fastp), val(organelle_type), path(db)
    )

    //
    // PHASE 3: Conditional reseed.
    //   Assess each first-pass result; for poor/failed assemblies (seed too
    //   divergent), download a closely related mitogenome and re-run seeded.
    //   Gated by params.skip_getorganelle_reseed (default: run).
    //

    // Join first-pass fasta + log + reads so we can assess and, if needed, reseed.
    ch_firstpass = GETORGANELLE_FROMREADS.out.fasta
        .join(GETORGANELLE_FROMREADS.out.log, by: 0)
        .join(fastp_with_mt_assembly_prefix, by: 0)   // [meta, fasta, log, reads]

    if (!params.skip_getorganelle_reseed) {

        ch_assessed = ch_firstpass.branch { _meta, fasta, log, _reads ->
            reseed: needsReseed(fasta, log)
            keep:   !needsReseed(fasta, log)
        }

        // Download a closely related mitogenome to use as the seed. Reuses
        // findMitoReference, now driven by the resolved reference_species_id.
        ch_reseed_meta = ch_assessed.reseed.map { meta, _fasta, _log, _reads -> meta }

        MITOHIFI_FINDMITOREFERENCE (
            ch_reseed_meta
        )

        ch_seed = MITOHIFI_FINDMITOREFERENCE.out.reference
            .map { meta, ref_fasta, _ref_gb -> [meta, ref_fasta] }   // [meta, seed]

        ch_reseed_reads = ch_assessed.reseed.map { meta, _fasta, _log, reads -> [meta, reads] }

        // Reseed candidates that actually got a seed downloaded.
        ch_reseed_input = ch_reseed_reads
            .join(ch_seed, by: 0)                       // [meta, reads, seed]
            .combine(ch_db)                             // [meta, reads, seed, org_type, db]
            .map { meta, reads, seed, org_type, db -> [meta, reads, org_type, db, seed] }

        GETORGANELLE_RESEED (
            ch_reseed_input // tuple val(meta), path(fastp), val(organelle_type), path(db), path(seed)
        )

        // Reseed candidates whose reference download was ignored/failed: fall
        // back to their (empty) first-pass result so no sample is dropped.
        ch_reseed_seedless = ch_assessed.reseed
            .map { meta, fasta, _log, _reads -> [meta, fasta] }
            .join(ch_seed, by: 0, remainder: true)
            .filter { it[2] == null }
            .map { meta, fasta, _seed -> [meta, fasta] }

        ch_reseed_seedless_log = ch_assessed.reseed
            .map { meta, _fasta, log, _reads -> [meta, log] }
            .join(ch_seed, by: 0, remainder: true)
            .filter { it[2] == null }
            .map { meta, log, _seed -> [meta, log] }

        // Final per-sample assembly = kept first-pass + reseed results + seedless fallback.
        ch_assembly_fasta = ch_assessed.keep.map { meta, fasta, _log, _reads -> [meta, fasta] }
            .mix(GETORGANELLE_RESEED.out.fasta)
            .mix(ch_reseed_seedless)
        ch_assembly_log = ch_assessed.keep.map { meta, _fasta, log, _reads -> [meta, log] }
            .mix(GETORGANELLE_RESEED.out.log)
            .mix(ch_reseed_seedless_log)

        // Reseed-specific files to fold into the collected channels below.
        ch_multiqc_files = ch_multiqc_files
            .mix(GETORGANELLE_RESEED.out.tool_params.collect { it[1] })
            .mix(MITOHIFI_FINDMITOREFERENCE.out.tool_params.collect { it[1] })
        ch_summary_files = ch_summary_files
            .mix(GETORGANELLE_RESEED.out.fasta.map { _meta, fasta -> fasta })
            .mix(GETORGANELLE_RESEED.out.log.map { _meta, log -> log })
            .mix(GETORGANELLE_RESEED.out.org_assm_graph)
            .mix(GETORGANELLE_RESEED.out.raw_assm_graph)
            .mix(GETORGANELLE_RESEED.out.simp_assm_graph)
        ch_versions = ch_versions
            .mix(GETORGANELLE_RESEED.out.versions.first())
            .mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())

    } else {
        ch_assembly_fasta = GETORGANELLE_FROMREADS.out.fasta
        ch_assembly_log   = GETORGANELLE_FROMREADS.out.log
    }

    //
    // Collect files
    //

    // ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_FROMREADS.out.etc.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_CONFIG.out.tool_params)
    ch_multiqc_files = ch_multiqc_files.mix(CAT_FASTQ.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_FROMREADS.out.tool_params.collect { it[1] })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.fasta.map { meta, fasta -> fasta })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.log.map { meta, log -> log })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.org_assm_graph)
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.raw_assm_graph)
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.simp_assm_graph)
    ch_versions = ch_versions.mix(GETORGANELLE_CONFIG.out.versions.first())
    ch_versions = ch_versions.mix(GETORGANELLE_FROMREADS.out.versions.first())
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // Emit outputs
    //

    emit:
    assembly_fasta  = ch_assembly_fasta
    assembly_log    = ch_assembly_log
    summary_files   = ch_summary_files
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
