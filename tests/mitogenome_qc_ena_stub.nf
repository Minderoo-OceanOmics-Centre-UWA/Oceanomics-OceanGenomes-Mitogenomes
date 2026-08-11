#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MITOGENOME_QC } from '../subworkflows/local/mitogenome_qc/main'

workflow {
    meta = [
        id: 'OG910',
        mt_assembly_prefix: 'OG910.hifi.250101.v3mitohifi',
        sequencing_type: 'hifi',
        date: '250101',
        genetic_code: 2
    ]
    annotations = file("${projectDir}/test_data/annotation_stub/*", checkIfExists: true)
    MITOGENOME_QC(
        channel.of(
            tuple(meta, 'Test species', true, true, annotations)
        )
    )
}
