process ENA_VALIDATION_RESULT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(validation_files, stageAs: "validation_inputs/*")
    val settings

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.ena_validation_result.tsv"), emit: record
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def webin_arg = settings.webin_requested ? '--webin-requested' : ''
    // Study is per candidate (one ENA child study per technology), so it comes
    // from meta rather than the run-level settings map.
    def ena_study = meta.ena_study?.toString()?.trim() ?: ''
    """
    collate_ena_validation.py record \\
        --input 'validation_inputs/*' \\
        --output '${meta.mt_assembly_prefix}.ena_validation_result.tsv' \\
        --assembly-prefix '${meta.mt_assembly_prefix}' \\
        --og-id '${meta.id}' \\
        --ena-study '${ena_study}' \\
        --validation-mode '${settings.validation_mode}' \\
        --validation-attempt '${settings.validation_attempt}' \\
        ${webin_arg}

    printf '"%s":\n    python: "%s"\n    collate_ena_validation: "1.0.0"\n' \\
        "${task.process}" "\$(python --version | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    printf 'assembly_prefix\tog_id\ttech\tseq_date\tcode\tena_study\tvalidation_mode\tvalidation_attempt\ttable2asn_status\treject_count\terror_count\twarning_count\tinfo_count\tfatal_discrepancy_count\tnostop_count\tblocking_codes\twarning_codes\tconversion_status\tconversion_reason\tconversion_exit\tpreflight_status\tpreflight_reason\tpreflight_exit\twebin_status\twebin_reason\twebin_exit\tsubmission_ready\n' > '${meta.mt_assembly_prefix}.ena_validation_result.tsv'
    printf '${meta.mt_assembly_prefix}\t${meta.id}\t\t\t\t${meta.ena_study ?: ''}\t${settings.validation_mode}\t${settings.validation_attempt}\tPASS\t0\t0\t0\t0\t0\t0\t\t\tPASS\tok\t0\tNOT_APPLICABLE\tnot_applicable\t\tPASS\tvalidated\t0\ttrue\n' >> '${meta.mt_assembly_prefix}.ena_validation_result.tsv'
    printf '"%s":\n    python: "stub"\n    collate_ena_validation: "stub"\n' "${task.process}" > versions.yml
    """
}
