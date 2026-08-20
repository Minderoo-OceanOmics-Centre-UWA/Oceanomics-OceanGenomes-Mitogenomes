nextflow.enable.dsl = 2

include { enaStudyAnnotate } from '../../subworkflows/local/utils_ena_targets/main'

workflow ENA_TARGETS {
    take:
    candidates

    main:
    resolved = candidates.map { meta -> enaStudyAnnotate(params, meta) }

    emit:
    resolved
}
