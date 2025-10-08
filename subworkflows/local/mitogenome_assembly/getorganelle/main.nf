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
include { PUSH_MTDNA_ASSM_RESULTS   } from '../../../../modules/local/upload_results/mtdna'

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
    //

    combined_input = fastp_with_mt_assembly_prefix.combine(GETORGANELLE_CONFIG.out.db)

    //
    // MODULE: Run assembly using GetOrganelle from reads
    //

    GETORGANELLE_FROMREADS (
        combined_input // tuple val(meta), path(fastp), val(organelle_type), path(db)
    )

    //
    // Collect files
    //

    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_FROMREADS.out.etc.collect{it[1]})
    ch_versions = ch_versions.mix(GETORGANELLE_CONFIG.out.versions.first())
    ch_versions = ch_versions.mix(GETORGANELLE_FROMREADS.out.versions.first())
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // Emit outputs
    //

    emit:
    assembly_fasta  = GETORGANELLE_FROMREADS.out.fasta
    assembly_log    = GETORGANELLE_FROMREADS.out.log
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
