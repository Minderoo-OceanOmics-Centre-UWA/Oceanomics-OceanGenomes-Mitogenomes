/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'
include { enaTargetAnnotate; validateEnaTargets } from '../utils_ena_targets/main'

// Mitogenome assembly
include { BUILD_SOURCE_MODIFIERS} from '../../../modules/local/genome_qc/build_source_modifiers'
include { FORMAT_FILES          } from '../../../modules/local/genome_qc/format_files'
include { ALLOCATE_ENA_LOCUS_TAGS } from '../../../modules/local/genome_qc/allocate_ena_locus_tags'
include { EXTRACT_GENES_GFF     } from '../../../modules/local/genome_qc/extract_genes/gff'
include { EXTRACT_GENES_GB      } from '../../../modules/local/genome_qc/extract_genes/gb'
include { TRANSLATE_GENES       } from '../../../modules/local/genome_qc/translate_genes'
include { GEN_FILES_TABLE2ASN   } from '../../../modules/local/genome_qc/gen_files_table2asn'
include { PARSE_TABLE2ASN_VALIDATION } from '../../../modules/local/genome_qc/parse_table2asn_validation'
include { ENA_FLATFILE          } from '../../../modules/local/genome_qc/ena_flatfile'
include { PREPARE_ENA_METADATA  } from '../../../modules/local/genome_qc/prepare_ena_metadata'
include { BUILD_ENA_CANDIDATE_PACKAGE } from '../../../modules/local/genome_qc/build_ena_candidate_package'
include { WEBIN_VALIDATE_GENOME } from '../../../modules/local/genome_qc/webin_validate_genome'
include { WEBIN_VALIDATE        } from '../../../modules/local/genome_qc/webin_validate'
include { ENA_VALIDATION_RESULT } from '../../../modules/local/genome_qc/ena_validation_result'
include { ENA_VALIDATION_SUMMARY} from '../../../modules/local/genome_qc/ena_validation_summary'
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

    // Fail on a mis-set study/prefix pair before any candidate is packaged, not
    // at the point one is about to be submitted under the wrong namespace.
    validateEnaTargets(params)

    // Per-sample circularity verdict (resolved at the QC gate) carried as a value
    // so the table2asn topology/completeness modifiers reflect the real assembly,
    // rather than asserting a circular topology on every sample. Dropped from the
    // FORMAT_FILES input so that process keeps its existing signature.
    // The four-part mt_assembly_prefix remains the output-directory identity.
    // The annotation-bearing FASTA/TBL basename is the public ENA sequence and
    // ASSEMBLYNAME identity (for example ...v3mitohifi.emma102).
    ch_ena_named = mitogenome_qc.map { meta, species, proceed, circular, files ->
        def annotation_tbl = files.find { it.name.endsWith('.tbl') }
        def annotation_fa = files.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') }
        def full_seqid = annotation_tbl?.baseName ?: annotation_fa?.baseName
        if (!full_seqid) {
            error "Cannot derive the annotated full SeqID for ${meta.id}: annotation bundle has no TBL/FASTA"
        }
        if (!full_seqid.startsWith("${meta.id}.")) {
            error "Annotated SeqID '${full_seqid}' does not belong to specimen ${meta.id}"
        }
        def annotation_version = full_seqid.tokenize('.').last()
        // The ENA study and its registered locus-tag prefix are properties of the
        // candidate's sequencing technology, resolved once here so every ENA step
        // downstream reads the same value from meta and no join key changes shape.
        def meta_ext = enaTargetAnnotate(params, meta + [
            full_seqid: full_seqid,
            annotation_version: annotation_version,
            scientific_name: species
        ])
        [ meta_ext, species, proceed, circular, files ]
    }

    ch_circular = ch_ena_named.map { meta, _species, _proceed, circular, _files -> [ meta, circular ] }
    ch_format_input = ch_ena_named.map { meta, species, proceed, _circular, files -> [ meta, species, proceed, files ] }

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

    // Allocate/reuse specimen-level gene serials before table2asn.  Serials are
    // shared across every technology for a specimen; the tag is rendered from
    // meta.ena_locus_prefix, so the same gene carries a different tag in each
    // technology's record, as ENA requires.  Every viable candidate version gets
    // a complete package; per-technology selection happens later over those.
    ALLOCATE_ENA_LOCUS_TAGS(
        FORMAT_FILES.out.gb,
        sql_config_file
    )
    ch_multiqc_files = ch_multiqc_files.mix(ALLOCATE_ENA_LOCUS_TAGS.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(ALLOCATE_ENA_LOCUS_TAGS.out.versions.first())

    PREPARE_ENA_METADATA(
        ch_ena_named.map { meta, _species, _proceed, _circular, _files -> meta },
        sql_config_file
    )
    ch_versions = ch_versions.mix(PREPARE_ENA_METADATA.out.versions.first())


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
        .join(ALLOCATE_ENA_LOCUS_TAGS.out.tagged_tbl, by: 0)
        .map { meta, sample_fa, _untagged_tbl, sample_cmt, tagged_tbl ->
            [ meta, sample_fa, tagged_tbl, sample_cmt ]
        }
        .join(BUILD_SOURCE_MODIFIERS.out.src_file, by:0)
        .join(ch_circular, by:0)

    GEN_FILES_TABLE2ASN (
        ch_processed_files, // tuple val(meta), path("processed/*.{fa,fasta}"), path("processed/*.{gb,tbl}"), path("processed/*.cmt"), path("*.src"), val(circular)
        template_sbt_file // sbt template generated from genbank, specific for OceanOmics
    )
    ch_table2asn_parser_input = GEN_FILES_TABLE2ASN.out.val_file
        .join(GEN_FILES_TABLE2ASN.out.discrepancy_file, by: 0)
        .join(ch_circular, by: 0)
    PARSE_TABLE2ASN_VALIDATION(ch_table2asn_parser_input)
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(GEN_FILES_TABLE2ASN.out.versions.first())
    // Feed raw and normalised table2asn validation output into MultiQC inputs.
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.val_file.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.stats_file.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GEN_FILES_TABLE2ASN.out.discrepancy_file.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PARSE_TABLE2ASN_VALIDATION.out.findings.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PARSE_TABLE2ASN_VALIDATION.out.status.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PARSE_TABLE2ASN_VALIDATION.out.qc_flags.collect { it[1] })
    ch_versions = ch_versions.mix(PARSE_TABLE2ASN_VALIDATION.out.versions.first())

    // ERROR/REJECT validator findings and FATAL discrepancy findings quarantine
    // only that sample. Warnings remain visible but continue to ENA conversion.
    ch_table2asn_pass = GEN_FILES_TABLE2ASN.out.gbf_file
        .join(PARSE_TABLE2ASN_VALIDATION.out.status, by: 0)
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

    ch_ena_package_input = FORMAT_FILES.out.processed_files
        .map { meta, sample_fa, _sample_tbl, _sample_cmt -> [ meta, sample_fa ] }
        // GFF only: FORMAT_FILES.out.gff also carries the fasta, and staging the
        // same basename twice would collide in the task work directory.
        .join(FORMAT_FILES.out.gff.map { meta, _sample_fa, sample_gff -> [ meta, sample_gff ] }, by: 0)
        .join(ENA_FLATFILE.out.embl_file, by: 0)
        .join(ALLOCATE_ENA_LOCUS_TAGS.out.mapping, by: 0)
        .join(ALLOCATE_ENA_LOCUS_TAGS.out.tagged_tbl, by: 0)
        .join(PREPARE_ENA_METADATA.out.metadata, by: 0)
    BUILD_ENA_CANDIDATE_PACKAGE(ch_ena_package_input)
    ch_multiqc_files = ch_multiqc_files
        .mix(BUILD_ENA_CANDIDATE_PACKAGE.out.validation.collect { it[1] })
        .mix(BUILD_ENA_CANDIDATE_PACKAGE.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(BUILD_ENA_CANDIDATE_PACKAGE.out.versions.first())

    ch_webin_test_status = Channel.empty()
    if (params.ena_validate_webin_test) {
        if (!secrets.WEBIN_USERNAME || !secrets.WEBIN_PASSWORD) {
            error "Nextflow secrets WEBIN_USERNAME and WEBIN_PASSWORD are required for ENA Webin test validation."
        }
        ch_ready_ena_packages = BUILD_ENA_CANDIDATE_PACKAGE.out.package_dir
            .join(BUILD_ENA_CANDIDATE_PACKAGE.out.metadata, by: 0)
            .filter { _meta, _package, metadata ->
                metadata.text.contains('"package_status": "READY"') &&
                    metadata.text.contains('"local_validation_status": "PASS"')
            }
            .map { meta, package_dir, _metadata -> [ meta, package_dir ] }
        WEBIN_VALIDATE_GENOME(
            ch_ready_ena_packages,
            'test',
            params.ena_validation_attempt
        )
        ch_webin_test_status = WEBIN_VALIDATE_GENOME.out.status
        ch_multiqc_files = ch_multiqc_files
            .mix(WEBIN_VALIDATE_GENOME.out.status.collect { it[2] })
            .mix(WEBIN_VALIDATE_GENOME.out.tool_params.collect { it[1] })
        ch_versions = ch_versions.mix(WEBIN_VALIDATE_GENOME.out.versions.first())
    }
    if (params.ena_validate_webin_production) {
        log.warn 'Production Webin validation runs on the selected package per technology in the standalone ENA selection workflow, not for every pipeline candidate.'
    }

    if (params.ena_webin_validate) {
        if (!secrets.WEBIN_USERNAME || !secrets.WEBIN_PASSWORD) {
            error "Nextflow secrets WEBIN_USERNAME and WEBIN_PASSWORD are required when --ena_webin_validate is enabled."
        }

        WEBIN_VALIDATE(ENA_FLATFILE.out.embl_file, params.ena_validation_attempt)
        ch_multiqc_files = ch_multiqc_files.mix(WEBIN_VALIDATE.out.tool_params.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(WEBIN_VALIDATE.out.status.collect { it[1] })
        ch_versions = ch_versions.mix(WEBIN_VALIDATE.out.versions.first())
        ch_validated_ena_flatfile = WEBIN_VALIDATE.out.validated_flatfile
    }

    // Collate every reached gate into one durable record per assembly. Missing
    // downstream files become explicit NOT_RUN/NOT_REQUESTED values rather than
    // silently dropping a quarantined sample.
    ch_ena_validation_files = PARSE_TABLE2ASN_VALIDATION.out.status
        .mix(ENA_FLATFILE.out.status)
        .mix(ENA_FLATFILE.out.checks)
        .mix(ENA_FLATFILE.out.embl_file)
        .mix(BUILD_ENA_CANDIDATE_PACKAGE.out.metadata)
        .mix(BUILD_ENA_CANDIDATE_PACKAGE.out.validation)
        .mix(ch_webin_test_status.map { meta, _service, status -> [meta, status] })
    if (params.ena_webin_validate) {
        ch_ena_validation_files = ch_ena_validation_files
            .mix(WEBIN_VALIDATE.out.status)
            .mix(WEBIN_VALIDATE.out.manifest)
    }

    // NOTE: a bare groupTuple, so no sample's record is built until every sample has cleared
    // QC. That is tolerable ONLY because this is the tail of the subworkflow -- table2asn,
    // the flatfile conversion, the package build and Webin all run per sample and are not
    // held by it, and ENA_VALIDATION_SUMMARY below collects across the whole run anyway.
    // Contrast the assembly-upload channel in the parent workflow, where the same operator
    // sat mid-pipeline and stopped every finished sample from reaching QC at all.
    //
    // Do NOT "fix" this with a groupKey carrying a constant size: the per-sample file count
    // genuinely varies. ENA_FLATFILE only runs on samples that passed table2asn
    // (ch_table2asn_pass above), BUILD_ENA_CANDIDATE_PACKAGE only on those whose whole join
    // chain matched, ch_webin_test_status only on READY+PASS packages, and WEBIN_VALIDATE
    // only under --ena_webin_validate. A too-large size makes groupTuple never emit for a
    // quarantined sample and silently lose it; a too-small one emits an incomplete record.
    // Sizing this correctly means first making each optional contributor emit an explicit
    // NOT_RUN row per sample -- which is the shape the comment above already assumes.
    ch_ena_validation_inputs = ch_ena_validation_files
        .map { meta, validation_file -> tuple(meta.mt_assembly_prefix, meta, validation_file) }
        .groupTuple(by: 0)
        .map { _prefix, metas, validation_files -> tuple(metas[0], validation_files.flatten()) }

    // No ena_study here: it is per candidate now and read from meta.ena_study.
    ena_validation_settings = [
        validation_mode: 'pipeline',
        validation_attempt: params.ena_validation_attempt,
        webin_requested: params.ena_webin_validate as boolean,
        workflow_run_name: workflow.runName ?: '',
        workflow_session_id: workflow.sessionId?.toString() ?: '',
        pipeline_revision: workflow.revision ?: workflow.commitId ?: ''
    ]
    ENA_VALIDATION_RESULT(ch_ena_validation_inputs, ena_validation_settings)
    ENA_VALIDATION_SUMMARY(ENA_VALIDATION_RESULT.out.record.map { _meta, record -> record }.collect())

    ch_multiqc_files = ch_multiqc_files.mix(ENA_VALIDATION_SUMMARY.out.multiqc)
    ch_versions = ch_versions
        .mix(ENA_VALIDATION_RESULT.out.versions.first())
        .mix(ENA_VALIDATION_SUMMARY.out.versions)
    
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
    ena_validation_records = ENA_VALIDATION_RESULT.out.record
    ena_validation_summary = ENA_VALIDATION_SUMMARY.out.multiqc
    ena_run_summary         = ENA_VALIDATION_SUMMARY.out.run_summary
    ena_candidate_packages  = BUILD_ENA_CANDIDATE_PACKAGE.out.package_dir
    ena_candidate_metadata  = BUILD_ENA_CANDIDATE_PACKAGE.out.metadata
    ena_webin_test_status   = ch_webin_test_status
}
