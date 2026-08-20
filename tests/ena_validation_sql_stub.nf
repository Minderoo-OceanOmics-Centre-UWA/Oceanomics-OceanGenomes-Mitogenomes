#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ENA_VALIDATION_RESULT  } from '../modules/local/genome_qc/ena_validation_result/main'
include { ENA_VALIDATION_SUMMARY } from '../modules/local/genome_qc/ena_validation_summary/main'
include { UPLOAD_ENA_RESULTS     } from '../subworkflows/local/upload_results_mito/main'

workflow {
    meta = [
        id: 'ena_stub',
        mt_assembly_prefix: 'ena_stub.hifi.260101.final',
        full_seqid: 'ena_stub.hifi.260101.final.emma102',
        ena_study: 'PRJEB123419'
    ]
    status = file("${projectDir}/test_data/ena_stub.table2asn_status.tsv", checkIfExists: true)
    flatfile = file("${projectDir}/test_data/ena_stub.embl.gz", checkIfExists: true)
    config = file("${projectDir}/../test_data/sql_config.txt", checkIfExists: true)
    inputs = channel.of(tuple(meta, [status, flatfile]))
    settings = [
        validation_mode: 'pipeline',
        validation_attempt: 'stub', webin_requested: false
    ]

    ENA_VALIDATION_RESULT(inputs, settings)
    ENA_VALIDATION_SUMMARY(ENA_VALIDATION_RESULT.out.record.map { _meta, record -> record }.collect())
    UPLOAD_ENA_RESULTS(ENA_VALIDATION_RESULT.out.record, Channel.empty(), config)
}
