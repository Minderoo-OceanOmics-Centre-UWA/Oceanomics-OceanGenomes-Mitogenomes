#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ENA_VALIDATION_RESULT  } from '../modules/local/genome_qc/ena_validation_result/main'
include { ENA_VALIDATION_SUMMARY } from '../modules/local/genome_qc/ena_validation_summary/main'
include { UPLOAD_ENA_RESULTS     } from '../subworkflows/local/upload_results_mito/main'

workflow {
    meta = [id: 'ena_stub', mt_assembly_prefix: 'ena_stub.hifi.260101.final']
    status = file("${projectDir}/test_data/ena_stub.table2asn_status.tsv", checkIfExists: true)
    flatfile = file("${projectDir}/test_data/ena_stub.embl.gz", checkIfExists: true)
    config = file("${projectDir}/../test_data/sql_config.txt", checkIfExists: true)
    inputs = channel.of(tuple(meta, [status, flatfile]))
    settings = [
        ena_study: 'PRJEB000000', validation_mode: 'pipeline',
        validation_attempt: 'stub', webin_requested: false,
        workflow_run_name: 'stub', workflow_session_id: 'stub', pipeline_revision: 'stub'
    ]

    ENA_VALIDATION_RESULT(inputs, settings)
    ENA_VALIDATION_SUMMARY(ENA_VALIDATION_RESULT.out.record.map { _meta, record -> record }.collect())
    UPLOAD_ENA_RESULTS(ENA_VALIDATION_RESULT.out.record, Channel.empty(), config)
}
