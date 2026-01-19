#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oceangenomesmitogenomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/oceangenomesmitogenomes
    Website: https://nf-co.re/oceangenomesmitogenomes
    Slack  : https://nfcore.slack.com/channels/oceangenomesmitogenomes
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { OCEANGENOMESMITOGENOMES } from './workflows/oceangenomesmitogenomes'
include { PREPARE_SAMPLESHEET     } from './subworkflows/local/prepare_samplesheet/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oceangenomesmitogenomes_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oceangenomesmitogenomes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_OCEANGENOMESMITOGENOMES {

    take:
    getorg_input    // tuple for getorganelle - tuple(meta, reads)
    mitohifi_input // tuple for mitohifi - tuple(meta, reads)
    
    main:

    // 
    // WORKFLOW: Run pipeline
    //
    OCEANGENOMESMITOGENOMES (
        getorg_input,    // tuple for getorganelle - tuple(meta, reads)
        mitohifi_input, // tuple for mitohifi - tuple(meta, reads)
        params.curated_blast_db,
        params.sql_config,
        params.organelle_type,
    )
    emit:
    multiqc_report = OCEANGENOMESMITOGENOMES.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    
    PREPARE_SAMPLESHEET (
        params.input,
        params.input_dir,
        params.samplesheet_prefix ?: "samplesheet"
    )
    
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        PREPARE_SAMPLESHEET.out.samplesheet
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_OCEANGENOMESMITOGENOMES (
        PREPARE_SAMPLESHEET.out.getorg_input,
        PREPARE_SAMPLESHEET.out.mitohifi_input
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_OCEANGENOMESMITOGENOMES.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
