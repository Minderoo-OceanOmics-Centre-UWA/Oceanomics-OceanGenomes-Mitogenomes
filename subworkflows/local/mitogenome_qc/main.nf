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
include { ENA_FLATFILE          } from '../../../modules/local/genome_qc/ena_flatfile'
include { WEBIN_VALIDATE        } from '../../../modules/local/genome_qc/webin_validate'
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
    mitogenome_qc // tuple val(meta), val(species_name), val(proceed_qc true/false), val(circular true/false), path(annotation/*)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_validated_ena_flatfile = channel.empty()
    def sql_config_file = file(params.sql_config, checkIfExists: true)
    def template_sbt_file = file(params.template_sbt, checkIfExists: true)

    // Per-sample circularity verdict (resolved at the QC gate) carried as a value
    // so the table2asn topology/completeness modifiers reflect the real assembly,
    // rather than asserting a circular topology on every sample. Dropped from the
    // FORMAT_FILES input so that process keeps its existing signature.
    ch_circular = mitogenome_qc.map { meta, _species, _proceed, circular, _files -> [ meta, circular ] }
    ch_format_input = mitogenome_qc.map { meta, species, proceed, _circular, files -> [ meta, species, proceed, files ] }

    //
    // MODULE: Format files to align better with GenBank requirements and generate the cmt file
    //
    // mitogenome_qc.view { "Input to FORMAT_FILES: $it" }

    FORMAT_FILES (
        ch_format_input   // tuple val(meta), val(species_name), val(proceed_qc true/false), path(input_dir) - input dir is the annotation outputs directory
    )
    // output is tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), path("processed/*.cmt"), emit: processed_files
    ch_multiqc_files = ch_multiqc_files.mix(FORMAT_FILES.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(FORMAT_FILES.out.versions.first())


    // meta_only_ch = mitogenome_qc.map { meta, species_name, proceed_qc, emma_path -> 
    //     meta 
    // }

    //
    // MODULE: Build the Source Modifiers table from SQL database
    //

    BUILD_SOURCE_MODIFIERS (
        FORMAT_FILES.out.meta, // val(meta)
        sql_config_file // val(db_config)
    )
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_SOURCE_MODIFIERS.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(BUILD_SOURCE_MODIFIERS.out.versions.first())


    //
    // MODULE: Extract all gene sequences including tRNAs, either using a gff file. 
    //

    EXTRACT_GENES_GFF (
        FORMAT_FILES.out.gff    //tuple val(meta), path(fasta), path(gff)
    )
    ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_GENES_GFF.out.tool_params.collect { it[1] })
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
    ch_multiqc_files = ch_multiqc_files.mix(TRANSLATE_GENES.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(TRANSLATE_GENES.out.versions.first())
    // ch_multiqc_files = ch_multiqc_files.mix(TRANSLATE_GENES.out.proteins_dir)
    
    //
    // MODULE: Generate files and run table2asn
    //
    ch_processed_files = FORMAT_FILES.out.processed_files
        .join(BUILD_SOURCE_MODIFIERS.out.src_file, by:0)
        .join(ch_circular, by:0)

    GEN_FILES_TABLE2ASN (
        ch_processed_files, // tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), path("processed/*.cmt"), path("*.src"), val(circular)
        template_sbt_file // sbt template generated from genbank, specific for OceanOmics
    )
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(GEN_FILES_TABLE2ASN.out.versions.first())
    // Feed raw and normalised table2asn validation output into MultiQC inputs.
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.val_file.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.stats_file.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.discrepancy_file.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.validation_findings.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.validation_status.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.qc_flags.collect { it[1] })

    // ERROR/REJECT validator findings and FATAL discrepancy findings quarantine
    // only that sample. Warnings remain visible but continue to ENA conversion.
    ch_table2asn_pass = GEN_FILES_TABLE2ASN.out.gbf_file
        .join(GEN_FILES_TABLE2ASN.out.validation_status, by: 0)
        .filter { _meta, _gbf, status_file ->
            def rows = status_file.readLines()
            rows.size() > 1 && rows[1].split('\\t', -1)[1] == 'PASS'
        }
        .map { meta, gbf, _status_file -> tuple(meta, gbf) }

    ENA_FLATFILE(ch_table2asn_pass)
    ch_multiqc_files = ch_multiqc_files.mix(ENA_FLATFILE.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ENA_FLATFILE.out.status.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ENA_FLATFILE.out.checks.collect { it[1] })
    ch_versions = ch_versions.mix(ENA_FLATFILE.out.versions.first())

    if (params.ena_webin_validate) {
        if (!params.ena_study) {
            error "--ena_study is required when --ena_webin_validate is enabled."
        }
        if (!secrets.WEBIN_USERNAME || !secrets.WEBIN_PASSWORD) {
            error "Nextflow secrets WEBIN_USERNAME and WEBIN_PASSWORD are required when --ena_webin_validate is enabled."
        }

        WEBIN_VALIDATE(ENA_FLATFILE.out.embl_file, params.ena_study, params.ena_validation_attempt)
        ch_multiqc_files = ch_multiqc_files.mix(WEBIN_VALIDATE.out.tool_params.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(WEBIN_VALIDATE.out.status.collect { it[1] })
        ch_versions = ch_versions.mix(WEBIN_VALIDATE.out.versions.first())
        ch_validated_ena_flatfile = WEBIN_VALIDATE.out.validated_flatfile
    }
    
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
    ena_flatfile            = ENA_FLATFILE.out.embl_file
    validated_ena_flatfile  = ch_validated_ena_flatfile
}
