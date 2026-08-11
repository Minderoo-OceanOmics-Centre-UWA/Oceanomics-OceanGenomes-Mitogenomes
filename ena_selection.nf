#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SELECT_ENA_SUBMISSION  } from './modules/local/genome_qc/select_ena_submission'
include { WEBIN_VALIDATE_GENOME } from './modules/local/genome_qc/webin_validate_genome'
include { RECORD_ENA_PACKAGE_VALIDATION } from './modules/local/genome_qc/record_ena_package_validation'

workflow {
    if (!params.ena_package_metadata) {
        error "--ena_package_metadata is required and must be a glob matching candidate *.package_metadata.json files."
    }
    if (!(params.ena_selection_mode in ['report', 'apply'])) {
        error "--ena_selection_mode must be report or apply."
    }

    ch_metadata = channel.fromPath(params.ena_package_metadata, checkIfExists: true)
        .collect()
        .ifEmpty { error "No ENA package metadata files matched: ${params.ena_package_metadata}" }

    decision_file = params.ena_decision_file ?
        file(params.ena_decision_file, checkIfExists: true) : []
    db_config = params.sql_config ?
        file(params.sql_config, checkIfExists: true) : []
    if (params.ena_selection_mode == 'apply' && !params.sql_config) {
        error "--sql_config is required when --ena_selection_mode apply."
    }

    SELECT_ENA_SUBMISSION(
        ch_metadata,
        params.ena_selection_mode,
        decision_file,
        db_config,
        params.ena_selected_by
    )

    if (params.ena_validate_webin_production) {
        if (params.ena_selection_mode != 'apply') {
            error "Production validation requires --ena_selection_mode apply so the selected package is registered first."
        }
        if (!secrets.WEBIN_USERNAME || !secrets.WEBIN_PASSWORD) {
            error "Nextflow secrets WEBIN_USERNAME and WEBIN_PASSWORD are required for production Webin validation."
        }
        ch_selected_packages = SELECT_ENA_SUBMISSION.out.selected
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row ->
                def package_dir = file(row.package_path, checkIfExists: true)
                tuple(
                    [
                        id: row.og_id,
                        full_seqid: row.full_seqid,
                        mt_assembly_prefix: row.full_seqid.tokenize('.')[0..3].join('.')
                    ],
                    package_dir
                )
            }
        WEBIN_VALIDATE_GENOME(
            ch_selected_packages,
            'production',
            params.ena_validation_attempt
        )
        RECORD_ENA_PACKAGE_VALIDATION(
            WEBIN_VALIDATE_GENOME.out.status,
            db_config
        )
    }
}
