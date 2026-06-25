/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data upload modules
include { PUSH_MTDNA_ASSM_RESULTS       } from '../../../modules/local/upload_results/mtdna'
include { SPECIES_VALIDATION            } from '../../../modules/local/species_validation'
include { PUSH_MTDNA_ANNOTATION_RESULTS } from '../../../modules/local/upload_results/emma'
include { PUSH_LCA_BLAST_RESULTS        } from '../../../modules/local/upload_results/lca'
include { PUSH_LCA_RAW_RESULTS          } from '../../../modules/local/upload_results/lca_raw'
include { UPLOAD_RESULTS_SUMMARY        } from '../../../modules/local/upload_results/summary'
include { EVALUATE_QC_CONDITIONS        } from '../../../modules/local/evaluate_qc_conditions'
include { QC_SUMMARY                    } from '../../../modules/local/multiqc/qc_summary'

// Helper functions
include { softwareVersionsToYAML        } from '../../nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DATE UPLOAD AND SPECIES CHECK WORKFLOW

        - Specific OceanOmics code using the PostgreSQL database.
        - Check the LCA results against nominal species ID and push results to SQL db
        - Determines species that have validated species ID to proceed with QC to prepare
          the sample for submittion to Genbank.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UPLOAD_RESULTS {

    take:
    assembly_results
    annotation_results
    blast_filtered_results
    lca_results
    lca_raw_results
    circularity_evidence // tuple val(meta), path(circularity_check.tsv) — HiFi only
    sql_config // params.sql_config

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 
    //
    // MODULE: Pulling out the statistics from the assembly files and updating the SQL database
    //

    PUSH_MTDNA_ASSM_RESULTS (
        assembly_results,
        params.sql_config
    )

    //
    // Re grouping the CO1, 12s and 16s BLAST and LCA results per sample
    //

    grouped_lca = lca_results
        .groupTuple(by: 0)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }

    grouped_blast = blast_filtered_results
        .groupTuple(by: 0)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }
    
    grouped_blast_lca = grouped_blast.join(grouped_lca, by: 0)
    
    //
    // MODULE: Checking the LCA results against the nominal species ID in the SQL database
    //

    SPECIES_VALIDATION (
        grouped_blast_lca, // tuple val(meta), path(blast_filtered), path(lca_filtered)
        sql_config // params.sql_config
    )

    //
    // MODULE: Calculating the statistics of the annotations and updating the SQL database
    //

    PUSH_MTDNA_ANNOTATION_RESULTS (
        annotation_results, // tuple val(meta), path("annotation/*")
        sql_config // params.sql_config
    )

    //
    // MODULE: Updating the SQL database with the LCA and filtered BLAST results
    //

    PUSH_LCA_BLAST_RESULTS (
        SPECIES_VALIDATION.out.full, // tuple path ("lca_combined.${meta.id}.tsv"), path ("blast_combined.${meta.id}.txt"),
        sql_config // params.sql_config
    )

    //
    // MODULE: Pushing raw per-hit LCA rows into the lca_raw_results table.
    //         Group the per-region lca_raw files by sample so each sample gets
    //         one upload call with all its regions.
    //
    grouped_lca_raw = lca_raw_results
        .groupTuple(by: 0)
        .map { tuple_ ->
            def meta = tuple_[0]
            def files = tuple_[1..-1].flatten()
            [meta, files]
        }

    PUSH_LCA_RAW_RESULTS (
        grouped_lca_raw, // tuple val(meta), path(lca_raw.*.tsv)
        sql_config // params.sql_config
    )

    //
    // SUBWORKFLOW: Evaluate the mitogenome and if it can continue on to final QC
    //
    /*  This solution provides a conditional QC subworkflow that evaluates two conditions and only proceeds with QC processes when **both** are satisfied:
        1. ✅ **LCA results**: Any row has the nominal species ID in the filtered blast
            results: `Found_in_blast_YN = "Yes"`
        2. ✅ **Annotation CSV**: The `passed` column value is `"yes"` when there are no
            missing genes and the genes are in the right order.
    */
    
    // Attach the per-sample circularity-check evidence (HiFi only) so the gate can
    // block length/repeat anomalies (concatemers, control-region VNTRs, unresolved
    // over-length assemblies) from progressing to QC. Samples without a check
    // (GetOrganelle, precomputed) get a "no anomaly" placeholder via a remainder
    // join, so they are never dropped and never blocked on this condition.
    def no_anomaly_file = file("${projectDir}/assets/empty_circularity_check.tsv", checkIfExists: true)
    ch_qc_conditions = SPECIES_VALIDATION.out.summary
        .join(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats, by: 0)
        .join(circularity_evidence, by: 0, remainder: true)
        // Keep only base rows (species + annotation present); drop evidence-only
        // remainder rows, which have a null in the blast slot.
        .filter { items -> items.size() >= 3 && items[1] != null }
        .map { items ->
            def circ = items.size() > 3 ? items[3] : null
            [items[0], items[1], items[2], circ ?: no_anomaly_file]
        }

    //
    // MODULE: evaluating the results to determine if to process the sample through QC
    //

    EVALUATE_QC_CONDITIONS (
        ch_qc_conditions // tuple val(meta), path(blast_filtered), path(annotation_stats.csv), path(circularity_check.tsv)
    )

    // Filter for samples that meet both conditions
    ch_qc_ready = EVALUATE_QC_CONDITIONS.out.evaluation
        .map { meta, species_file, proceed_file, circular_file ->
            def species_name = species_file.text.trim()
            def proceed_qc = proceed_file.text.trim()
            def circular = circular_file.text.trim()
            return [ meta, species_name, proceed_qc, circular ]
        }
        .filter { meta, species_name, proceed_qc, circular ->
            proceed_qc == "true"
        }

    // Log samples that will proceed to QC
    ch_qc_ready.view { meta, species_name, proceed_qc, circular ->
        "Sample ${meta.id} will proceed to QC with species: ${species_name} (circular: ${circular})"
    }

    // Filter for samples that dont meet the conditions
    ch_not_qc_ready = EVALUATE_QC_CONDITIONS.out.evaluation
        .map { meta, species_file, proceed_file, circular_file ->
            def species_name = species_file.text.trim()
            def proceed_qc = proceed_file.text.trim()
            def circular = circular_file.text.trim()
            return [ meta, species_name, proceed_qc, circular ]
        }
        .filter { meta, species_name, proceed_qc, circular ->
            proceed_qc == "false"
        }
        .view { meta, species_name, proceed_qc, circular ->
            "Sample ${meta.id} will NOT proceed to QC - conditions not met"
        }

    
    
    //
    // Build a per-sample QC summary TSV for MultiQC
    qc_summary_input = EVALUATE_QC_CONDITIONS.out.evaluation
        .join(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats, by: 0)
        .map { meta, species_file, proceed_file, circular_file, annotation_csv -> [ meta, species_file, proceed_file, annotation_csv ] }
    QC_SUMMARY (
        qc_summary_input // tuple val(meta), path(species_name.txt), path(proceed_qc.txt), path(annotation_stats.csv)
    )

    
    //
    // MODULE: Consolidate the per-step SQL upload status files into a single
    //         summary TSV + detailed appendix so all upload outcomes can be
    //         reviewed in one place instead of grep'ing many small files.
    //

    ch_upload_status_files = Channel.empty()
        .mix(PUSH_MTDNA_ASSM_RESULTS.out.upload)
        .mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.upload.map { _meta, f -> f })
        .mix(SPECIES_VALIDATION.out.upload.map { _meta, f -> f })
        .mix(PUSH_LCA_BLAST_RESULTS.out.upload)
        .mix(PUSH_LCA_RAW_RESULTS.out.upload)

    UPLOAD_RESULTS_SUMMARY (
        ch_upload_status_files.collect()
    )

    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    //  - Species validation outputs (per-sample TSVs)
    //  - Annotation statistics CSVs
    //  - QC evaluation flags (for quick visibility in report) and summary table
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_VALIDATION.out.summary.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.evaluation.map { meta, species_file, proceed_file, circular_file -> proceed_file })
    ch_multiqc_files = ch_multiqc_files.mix(QC_SUMMARY.out.table.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ASSM_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_VALIDATION.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_LCA_BLAST_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_LCA_RAW_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_RESULTS_SUMMARY.out.summary)
    ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_RESULTS_SUMMARY.out.tool_params)
    ch_versions = ch_versions.mix(PUSH_MTDNA_ASSM_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(SPECIES_VALIDATION.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_BLAST_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_RAW_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(EVALUATE_QC_CONDITIONS.out.versions)
    ch_versions = ch_versions.mix(UPLOAD_RESULTS_SUMMARY.out.versions)



    //
    // Emit outputs
    //

    emit:
    qc_ready    = ch_qc_ready                   // channel: [ val(meta), val(species_name), val(proceed_qc true/false), val(circular true/false) ]
    assembly_summary_files = PUSH_MTDNA_ANNOTATION_RESULTS.out.stats.map { meta, stats -> stats }
    upload_summary = UPLOAD_RESULTS_SUMMARY.out.summary    // channel: path(upload_results_summary.tsv)
    upload_appendix = UPLOAD_RESULTS_SUMMARY.out.appendix  // channel: path(upload_results_appendix.txt)
    multiqc_files = ch_multiqc_files            // channel: [ path(multiqc_files) ]
    versions = ch_versions             // channel: [ path(versions.yml) ]
}
