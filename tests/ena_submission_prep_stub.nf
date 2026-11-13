#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENA_SUBMISSION_PREP } from '../subworkflows/local/ena_submission_prep/main'

// ENA_SUBMISSION_PREP resolves placeholders from ${projectDir}/assets, and nextflow
// sets projectDir to the directory of the launched script -- this one. tests/assets is
// a symlink back to the repo's own, so the subworkflow loads exactly the placeholders
// it would in a real run. Same arrangement as tests/local_qc_no_db/assets.
workflow {
    meta = [
        id: 'OG910',
        mt_assembly_prefix: 'OG910.hifi.250101.v3mitohifi',
        sequencing_type: 'hifi',
        date: '250101',
        genetic_code: 2
    ]
    annotations = file("${projectDir}/test_data/annotation_stub/*", checkIfExists: true)
    // Header-only depth placeholder: the stub asserts wiring, not coverage values.
    no_depth = file("${projectDir}/../assets/placeholders/empty_mito_depth.tsv", checkIfExists: true)
    ENA_SUBMISSION_PREP(
        channel.of(
            tuple(meta, 'Test species', true, true, annotations, no_depth)
        )
    )
}
