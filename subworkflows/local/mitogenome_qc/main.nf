/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'

// Mitogenome assembly
include { BUILD_SOURCE_MODIFIERS} from '../../../modules/local/genome_qc/build_source_modifiers'
include { FORMAT_FILES          } from '../../../modules/local/genome_qc/format_files'
include { EXTRACT_GENES_GFF     } from '../../../modules/local/genome_qc/extract_genes/gff'
include { EXTRACT_GENES_GB      } from '../../../modules/local/genome_qc/extract_genes/gb'
include { TRANSLATE_GENES       } from '../../../modules/local/genome_qc/translate_genes'
include { GEN_FILES_TABLE2ASN   } from '../../../modules/local/genome_qc/gen_files_table2asn'
// include { DIAGNOSTICS           } from '../../../modules/local/genome_qc/diagnostics'
// include { GROUPER               } from '../../../modules/local/genome_qc/grouper'
// include { SUBMITTER             } from '../../../modules/local/genome_qc/submitter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ASSEMBLY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_QC {

    take:
    mitogenome_qc // tuple val(meta), val(species_name), val(proceed_qc true/false), path(emma/*)
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    //
    // MODULE: Format files to align better with GenBank requirements and generate the cmt file
    //
    mitogenome_qc.view { "Input to FORMAT_FILES: $it" }

    FORMAT_FILES (
        mitogenome_qc     // tuple val(meta), val(species_name), val(proceed_qc true/false), path(input_dir) - input dir is the annotation outputs directory
    )
    // output is tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), path("processed/*.cmt"), emit: processed_files
    ch_versions = ch_versions.mix(FORMAT_FILES.out.versions.first())


    // meta_only_ch = mitogenome_qc.map { meta, species_name, proceed_qc, emma_path -> 
    //     meta 
    // }

    //
    // MODULE: Build the Source Modifiers table from SQL database
    //

    BUILD_SOURCE_MODIFIERS (
        FORMAT_FILES.out.meta, // val(meta)
        params.sql_config // val(db_config)
    )
    ch_versions = ch_versions.mix(BUILD_SOURCE_MODIFIERS.out.versions.first())


    //
    // MODULE: Extract all gene sequences including tRNAs, either using a gff file. 
    //

    EXTRACT_GENES_GFF (
        FORMAT_FILES.out.gff    //tuple val(meta), path(fasta), path(gff)
    )
    ch_versions = ch_versions.mix(EXTRACT_GENES_GFF.out.versions.first())

    //
    // MODULE: Extract all coding sequenses using a tbl/gb file. 
    //

    // EXTRACT_GENES_GB (
    //     FORMAT_FILES.out.gb    //tuple val(meta), path(fasta), path(tbl)
    // )
    // ch_versions = ch_versions.mix(EXTRACT_GENES_GB.out.versions.first())

    //
    // MODULE: Translate all cds to protein sequence
    //

    TRANSLATE_GENES (
        EXTRACT_GENES_GFF.out.genes_dir    // tuple val(meta), path(genes_path)
    )
    ch_versions = ch_versions.mix(TRANSLATE_GENES.out.versions.first())
    // ch_multiqc_files = ch_multiqc_files.mix(TRANSLATE_GENES.out.proteins_dir)
    
    //
    // MODULE: Generate files and run table2asn
    //
    ch_processed_files = FORMAT_FILES.out.processed_files.join(BUILD_SOURCE_MODIFIERS.out.src_file, by:0)

    GEN_FILES_TABLE2ASN (
        ch_processed_files, // tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), path("processed/*.cmt"), path("*.src") 
        params.template_sbt // sbt template generated from genbank, specific for OceanOmics
    )
    ch_versions = ch_versions.mix(GEN_FILES_TABLE2ASN.out.versions.first())
    // Feed table2asn validation output into MultiQC inputs (text summary)
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.val_file.collect { it[1] })
    
    //
    // MODULE: Interperate the translation diagnostics from table2asn output
    //

    // DIAGNOSTICS (

    // )

    //
    // MODULE: Group similar mitogenomes
    //
    /* Need to write a module for here that will check the database for other mitogenomes
        that have not been submitted or have an accession number in the SQL database.
        It will then find other mitogenomes that have the same transflated protein
        fingerprint and group them for submission.
        Maybe the mitogenomes that are ready for submission need to be noted in the SQL
        database so that they can be grouped together. Maybe there can be a ready directory
        and mitogenomes are put into a directory with other ones that match them and at then
        periodically they get sent, either when they hit 10 or at the start of each week. */


    //
    // MODULE: Genbank submitter
    //

    // SUBMITTER (

    // )


    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    // Need to update this section to include everything
    // ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect{it[1]})
    // ch_versions = ch_versions.mix(EMMA.out.versions.first())
    // ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    // ch_versions = ch_versions.mix(LCA.out.versions.first())



    //
    // Emit outputs
    //

    emit:
    multiqc_files           = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions                = ch_versions              // channel: [ path(versions.yml) ]
}
