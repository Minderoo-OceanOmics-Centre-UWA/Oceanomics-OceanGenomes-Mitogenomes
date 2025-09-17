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
        .groupTuple(by: 0, size: 3)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }

    grouped_blast = blast_filtered_results
        .groupTuple(by: 0, size: 3)
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
        annotation_results, // tuple val(meta), path("emma/*")
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
    // SUBWORKFLOW: Evaluate the mitogenome and if it can continue on to final QC
    //
    /*  This solution provides a conditional QC subworkflow that evaluates two conditions and only proceeds with QC processes when **both** are satisfied:
        1. ✅ **LCA results**: Any row has the nominal species ID in the filtered blast
            results: `Found_in_blast_YN = "Yes"`
        2. ✅ **Annotation CSV**: The `passed` column value is `"yes"` when there are no
            missing genes and the genes are in the right order.
    */
    
    ch_qc_conditions = SPECIES_VALIDATION.out.summary.join(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats, by:0)

    //
    // MODULE: evaluating the results to determine if to process the sample through QC
    //

    EVALUATE_QC_CONDITIONS (
        ch_qc_conditions // tuple val(meta), path(lca_results.${meta.id}.tsv), path(${meta.id}.annotation_stats.csv)
    )

    // Filter for samples that meet both conditions
    ch_qc_ready = EVALUATE_QC_CONDITIONS.out.evaluation
        .map { meta, species_file, proceed_file ->
            def species_name = species_file.text.trim()
            def proceed_qc = proceed_file.text.trim()
            return [ meta, species_name, proceed_qc ]
        }
        .filter { meta, species_name, proceed_qc ->
            proceed_qc == "true"
        }
    
    // Log samples that will proceed to QC
    ch_qc_ready.view { meta, species_name, proceed_qc ->
        "Sample ${meta.id} will proceed to QC with species: ${species_name}"
    }
    
    // Filter for samples that dont meet the conditions
    ch_not_qc_ready = EVALUATE_QC_CONDITIONS.out.evaluation
        .map { meta, species_file, proceed_file ->
            def species_name = species_file.text.trim()
            def proceed_qc = proceed_file.text.trim()
            return [ meta, species_name, proceed_qc ]
        }
        .filter { meta, species_name, proceed_qc ->
            proceed_qc == "false"
        }
        .view { meta, species_name, proceed_qc ->
            "Sample ${meta.id} will NOT proceed to QC - conditions not met"
        }

    
    
    //
    // Build a per-sample QC summary TSV for MultiQC
    qc_summary_input = EVALUATE_QC_CONDITIONS.out.evaluation
        .join(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats, by: 0)
        .map { meta, species_file, proceed_file, annotation_csv -> [ meta, species_file, proceed_file, annotation_csv ] }
    QC_SUMMARY (
        qc_summary_input // tuple val(meta), path(species_name.txt), path(proceed_qc.txt), path(annotation_stats.csv)
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
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.evaluation.map { meta, species_file, proceed_file -> proceed_file })
    ch_multiqc_files = ch_multiqc_files.mix(QC_SUMMARY.out.table.collect { it[1] })
    ch_versions = ch_versions.mix(PUSH_MTDNA_ASSM_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(SPECIES_VALIDATION.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_BLAST_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(EVALUATE_QC_CONDITIONS.out.versions)



    //
    // Emit outputs
    //

    emit:
    qc_ready    = ch_qc_ready                   // channel: [ val(meta), val(species_name), val(proceed_qc true/false) ]
    multiqc_files = ch_multiqc_files            // channel: [ path(multiqc_files) ]
    versions = ch_versions             // channel: [ path(versions.yml) ]
}
