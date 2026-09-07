/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'
include { enaStudyAnnotate; validateEnaStudy } from '../utils_ena_targets/main'

// Mitogenome assembly
include { BUILD_SOURCE_MODIFIERS} from '../../../modules/local/genome_qc/build_source_modifiers'
include { FORMAT_FILES          } from '../../../modules/local/genome_qc/format_files'
include { EXTRACT_GENES_GFF     } from '../../../modules/local/genome_qc/extract_genes/gff'
include { EXTRACT_GENES_GB      } from '../../../modules/local/genome_qc/extract_genes/gb'
include { TRANSLATE_GENES       } from '../../../modules/local/genome_qc/translate_genes'
include { GEN_FILES_TABLE2ASN   } from '../../../modules/local/genome_qc/gen_files_table2asn'
include { PARSE_TABLE2ASN_VALIDATION } from '../../../modules/local/genome_qc/parse_table2asn_validation'
include { ENA_FLATFILE          } from '../../../modules/local/genome_qc/ena_flatfile'
include { PREPARE_ENA_METADATA  } from '../../../modules/local/genome_qc/prepare_ena_metadata'
include { BUILD_ENA_CANDIDATE_PACKAGE } from '../../../modules/local/genome_qc/build_ena_candidate_package'
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

    // Fail on a missing or malformed study before any candidate is packaged, not
    // at the point one is about to be submitted under the wrong namespace.
    validateEnaStudy(params)

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
        // The run's ENA study is attached once here so every ENA step downstream
        // reads the same value from meta and no join key changes shape.
        def meta_ext = enaStudyAnnotate(params, meta + [
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

    // INSDC locus tags are deliberately absent from everything this pipeline
    // emits: a separate submission pipeline owns tag allocation and injects them
    // into the flat file before Webin.  table2asn will report NO_LOCUS_TAGS for
    // every record as a result, which is why that discrepancy code is advisory
    // in bin/parse_table2asn_validation.py.
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
    // Branch (not filter) so the quarantined complement is available locally: it is the exact
    // set of samples ENA_FLATFILE / package / Webin never run for, and it lets each of those
    // optional stages emit an explicit per-sample NOT_RUN row below, making the ENA validation
    // record a FIXED-size group that releases per sample (no close-time remainder). gbf_file is
    // a required table2asn output, so pass + quarantined partition every sample exactly.
    ch_table2asn_branched = GEN_FILES_TABLE2ASN.out.gbf_file
        .join(PARSE_TABLE2ASN_VALIDATION.out.status, by: 0)
        .branch { _meta, _gbf, status_file ->
            def rows = status_file.readLines()
            pass: rows.size() > 1 && rows[1].split('\\t', -1)[1] == 'PASS'
            fail: true
        }
    ch_table2asn_pass = ch_table2asn_branched.pass.map { meta, gbf, _status_file -> tuple(meta, gbf) }
    // Quarantined metas: no ENA_FLATFILE, package or Webin runs for these.
    ch_ena_quarantined = ch_table2asn_branched.fail.map { meta, _gbf, _status_file -> meta }

    // Headerless per-sample fragments for the run-level held_samples.tsv: the
    // table2asn quarantine set, with the blocking validator codes (status.tsv
    // column 9). table2asn FAIL is terminal -- no feedback loop -- so surfacing
    // it here is the only record outside the ENA validation summary.
    ch_held_fragments = ch_table2asn_branched.fail
        .collectFile { meta, _gbf, status_file ->
            def lines = status_file.readLines()
            def cols = lines.size() > 1 ? lines[1].split('\t', -1) : []
            def blocking = (cols.size() > 8 && cols[8]?.trim()) ? cols[8].trim() : 'unknown'
            [ "${meta.mt_assembly_prefix}.TABLE2ASN.held.tsv",
              "${meta.id}\t${meta.mt_assembly_prefix}\tTABLE2ASN\tFAIL_TABLE2ASN: ${blocking}\n" ]
        }

    ENA_FLATFILE(ch_table2asn_pass)
    ch_multiqc_files = ch_multiqc_files.mix(ENA_FLATFILE.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ENA_FLATFILE.out.status.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ENA_FLATFILE.out.checks.collect { it[1] })
    ch_versions = ch_versions.mix(ENA_FLATFILE.out.versions.first())

    // ENA_FLATFILE writes the .embl(.gz) -- and therefore package + Webin run -- ONLY when
    // conversion_status == PASS (the module gzips the flat file solely on PASS). Branch its
    // always-emitted status so "converted" vs "conversion failed" is a local per-sample value,
    // giving the exact complement for the embl / package / Webin NOT_RUN rows below without a
    // remainder join. Samples that reached ENA_FLATFILE but did not convert have no embl,
    // package or Webin outputs, exactly like the quarantined ones.
    ch_flatfile_status_branched = ENA_FLATFILE.out.status.branch { _meta, status_file ->
        def rows = status_file.readLines()
        pass: rows.size() > 1 && rows[1].split('\\t', -1)[1] == 'PASS'
        fail: true
    }
    ch_ena_no_embl = ch_ena_quarantined
        .mix( ch_flatfile_status_branched.fail.map { meta, _status_file -> meta } )

    // Flat file format check first: it is the pipeline's last ENA gate, and the
    // package records its verdict, so the build has to see the result. The
    // status file is emitted for a failing flat file too, so a FAIL still
    // produces a package that says why -- it does not silently drop a sample.
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

    ch_ena_package_base = FORMAT_FILES.out.processed_files
        .map { meta, sample_fa, sample_tbl, _sample_cmt -> [ meta, sample_fa, sample_tbl ] }
        // GFF only: FORMAT_FILES.out.gff also carries the fasta, and staging the
        // same basename twice would collide in the task work directory.
        .join(FORMAT_FILES.out.gff.map { meta, _sample_fa, sample_gff -> [ meta, sample_gff ] }, by: 0)
        // Concatenated gene FASTA only: genes_dir would also stage the per-CDS singles.
        .join(EXTRACT_GENES_GFF.out.genes_fa, by: 0)
        .join(ENA_FLATFILE.out.embl_file, by: 0)
        .join(PREPARE_ENA_METADATA.out.metadata, by: 0)

    // An empty list stages no file, which the module reads as NOT_REQUESTED.
    // Without this the join would starve the build whenever validation is off.
    ch_ena_package_input = params.ena_webin_validate
        ? ch_ena_package_base.join(WEBIN_VALIDATE.out.status, by: 0)
        : ch_ena_package_base.map { it + [ [] ] }
    BUILD_ENA_CANDIDATE_PACKAGE(ch_ena_package_input)
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_ENA_CANDIDATE_PACKAGE.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(BUILD_ENA_CANDIDATE_PACKAGE.out.versions.first())

    // Collate every reached gate into one durable record per assembly. Each optional
    // contributor is made TOTAL by emitting an explicit per-sample NOT_RUN row for the samples
    // it did not run on, so EVERY sample contributes the SAME fixed number of files:
    //   PARSE status (always) + conversion status + conversion checks + embl + package metadata
    //   [+ Webin status + Webin manifest]  =  5  (or 7 with --ena_webin_validate).
    // With the count fixed and known, a constant groupKey releases each sample's record the
    // moment its own files arrive -- no bare groupTuple, and NO close-time remainder -- so a
    // quarantined or conversion-failed sample streams exactly like a fully-passing one.
    //
    // Complements are local per-sample sets (no remainder join): conversion status + checks run
    // for every table2asn PASS, so their NOT_RUN complement is the quarantined set; embl,
    // package and Webin run only when conversion PASSed, so their complement is quarantined +
    // conversion-failed (ch_ena_no_embl).
    //
    // The NOT_RUN placeholders are named so collate_ena_validation.py IGNORES them: it reads
    // only the .table2asn_status.tsv / .ena_conversion_status.tsv / .webin_status.tsv suffixes
    // and already infers NOT_RUN / SKIPPED / NOT_REQUESTED from a stage's ABSENCE. The record is
    // therefore byte-identical to the bare-groupTuple version -- collate sees the same real
    // status files and the same absences; the placeholders only make the group size fixed.
    def notrun_flatfile_status  = file("${projectDir}/assets/placeholders/ena_not_run/flatfile_status.not_run",  checkIfExists: true)
    def notrun_flatfile_checks  = file("${projectDir}/assets/placeholders/ena_not_run/flatfile_checks.not_run",  checkIfExists: true)
    def notrun_flatfile_embl    = file("${projectDir}/assets/placeholders/ena_not_run/flatfile_embl.not_run",    checkIfExists: true)
    def notrun_package_metadata = file("${projectDir}/assets/placeholders/ena_not_run/package_metadata.not_run", checkIfExists: true)
    def notrun_webin_status     = file("${projectDir}/assets/placeholders/ena_not_run/webin_status.not_run",     checkIfExists: true)
    def notrun_webin_manifest   = file("${projectDir}/assets/placeholders/ena_not_run/webin_manifest.not_run",   checkIfExists: true)

    ch_ena_validation_files = PARSE_TABLE2ASN_VALIDATION.out.status                              // always, all samples
        .mix( ENA_FLATFILE.out.status )                                                          // conversion status: PASS set
        .mix( ch_ena_quarantined.map { meta -> [ meta, notrun_flatfile_status ] } )              //   + NOT_RUN: quarantined
        .mix( ENA_FLATFILE.out.checks )                                                          // conversion checks: PASS set
        .mix( ch_ena_quarantined.map { meta -> [ meta, notrun_flatfile_checks ] } )              //   + NOT_RUN: quarantined
        .mix( ENA_FLATFILE.out.embl_file )                                                       // embl: converted set
        .mix( ch_ena_no_embl.map { meta -> [ meta, notrun_flatfile_embl ] } )                    //   + NOT_RUN: quarantined + conv-failed
        .mix( BUILD_ENA_CANDIDATE_PACKAGE.out.metadata )                                         // package metadata: converted set
        .mix( ch_ena_no_embl.map { meta -> [ meta, notrun_package_metadata ] } )                 //   + NOT_RUN: quarantined + conv-failed
    def ena_record_slots = (params.ena_webin_validate as boolean) ? 7 : 5
    if (params.ena_webin_validate) {
        ch_ena_validation_files = ch_ena_validation_files
            .mix( WEBIN_VALIDATE.out.status )                                                     // Webin status: converted set
            .mix( ch_ena_no_embl.map { meta -> [ meta, notrun_webin_status ] } )                  //   + NOT_RUN: quarantined + conv-failed
            .mix( WEBIN_VALIDATE.out.manifest )                                                   // Webin manifest: converted set
            .mix( ch_ena_no_embl.map { meta -> [ meta, notrun_webin_manifest ] } )                //   + NOT_RUN: quarantined + conv-failed
    }

    // Fixed, known contributor count per sample -> constant groupKey + plain groupTuple, so
    // every sample (passing, quarantined or conversion-failed) releases its record on its own
    // as soon as its `ena_record_slots` files arrive. ENA_VALIDATION_SUMMARY below still
    // collects across the whole run.
    ch_ena_validation_inputs = ch_ena_validation_files
        .map { meta, validation_file -> tuple( groupKey(meta.full_seqid ?: meta.mt_assembly_prefix, ena_record_slots), meta, validation_file ) }
        .groupTuple(by: 0)
        .map { _key, metas, validation_files -> tuple(metas[0], validation_files.flatten()) }

    // No ena_study here: it is per candidate now and read from meta.ena_study.
    ena_validation_settings = [
        validation_mode: 'pipeline',
        validation_attempt: params.ena_validation_attempt,
        webin_requested: params.ena_webin_validate as boolean
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
    held_fragments          = ch_held_fragments            // channel: path(<prefix>.held.tsv) — one row per table2asn quarantine
    ena_candidate_packages  = BUILD_ENA_CANDIDATE_PACKAGE.out.package_dir
    ena_candidate_metadata  = BUILD_ENA_CANDIDATE_PACKAGE.out.metadata
}
