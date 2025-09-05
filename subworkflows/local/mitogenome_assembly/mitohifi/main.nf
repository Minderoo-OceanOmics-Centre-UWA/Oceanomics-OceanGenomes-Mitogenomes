/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../../nf-core/utils_nfcore_pipeline'

// Mitogenome assembly
include { SPECIES_QUERY                         } from '../../../../modules/local/species_query'
include { MITOHIFI_FINDMITOREFERENCE       } from '../../../../modules/nf-core/mitohifi/findmitoreference'
include { CAT_FASTQ                         } from '../../../../modules/nf-core/cat/fastq'
include { MITOHIFI_MITOHIFI    } from '../../../../modules/nf-core/mitohifi/mitohifi'
include { PUSH_MTDNA_ASSM_RESULTS   } from '../../../../modules/local/upload_results/mtdna'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ASSEMBLY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ASSEMBLY_MITOHIFI {

    take:
    fastp_reads // tuple val(meta), path(fastp)
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // map just the meta for the species query
    //

    ch_species_query = fastp_reads
    .map { meta, files -> meta }
    .view()
    
    //
    // MODULE: Extract nominal species from database
    //

    SPECIES_QUERY (
        ch_species_query,
        params.sql_config
    )

    //
    // MODULE: Find a closely related species for reference
    //

    MITOHIFI_FINDMITOREFERENCE (
        SPECIES_QUERY.out.species
    )

    //
    // MODULE: Concatenate fastq reads where there are multiple fastq files
    //

    CAT_FASTQ (
        fastp_reads
    )

    //
    // Combine fastp files with the mito reference output
    //

    combined_ch = CAT_FASTQ.out.reads.join(MITOHIFI_FINDMITOREFERENCE.out.reference, by: 0)

    //
    // Get the MitoHifi version to be used in the naming. Ensures the proper version is always used.
    //

    version_ch = MITOHIFI_FINDMITOREFERENCE.out.versions_tuple
    .map { meta, versions_file ->
        def content = versions_file.text
        def pattern = /mitohifi:\s*([^\s\n\r]+)/
        def matcher = content =~ pattern
        def version = matcher ? matcher[0][1].replaceAll(/["']/, '') : "unknown"
        [meta, version]
    }

    //
    // Embed the assembly prefix into meta.
    //

    combined_with_mt_assembly_prefix = combined_ch.join(version_ch, by: 0)
    .map { meta, fasta, ref_fasta, ref_gb, version ->  // Destructure all 3 elements correctly
        def version_stripped = version.replaceAll('\\.', '')
        def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${version_stripped}mitohifi"
        def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
        [meta_ext, fasta, ref_fasta, ref_gb]  // Return only what you want
    }


    //
    // MODULE: Run assembly using MitoHifi from reads
    //

    MITOHIFI_MITOHIFI (
        combined_with_mt_assembly_prefix,
        "r",
        "2"  // mito_code - Organism genetic code following NCBI table (for mitogenome annotation) using 2. Vertebrate mitogenome
    )

    //
    // Collect files
    //

    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    ch_versions = ch_versions.mix(MITOHIFI_MITOHIFI.out.versions.first())


    //
    // Collate and save software versions
    //

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oceangenomes_draftgenomes_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // Emit outputs
    //

    emit:
    assembly_fasta  = MITOHIFI_MITOHIFI.out.fasta
    assembly_log    = MITOHIFI_MITOHIFI.out.logs
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_collated_versions              // channel: [ path(versions.yml) ]
}
