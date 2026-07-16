#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENA_EMBL_PREFLIGHT } from '../modules/local/genome_qc/ena_embl_preflight'

workflow {
    meta = [id: 'ena_stub', mt_assembly_prefix: 'ena_stub']
    embl = file("${projectDir}/test_data/ena_stub.embl.gz", checkIfExists: true)
    ENA_EMBL_PREFLIGHT(channel.of(tuple(meta, embl)))
}
