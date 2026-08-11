nextflow.enable.dsl = 2

include { enaTargetAnnotate } from '../../subworkflows/local/utils_ena_targets/main'

workflow ENA_TARGETS {
    take:
    candidates

    main:
    resolved = candidates.map { meta -> enaTargetAnnotate(params, meta) }

    emit:
    resolved
}
