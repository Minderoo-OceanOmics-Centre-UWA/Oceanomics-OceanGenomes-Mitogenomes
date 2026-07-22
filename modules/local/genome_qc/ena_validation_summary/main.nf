process ENA_VALIDATION_SUMMARY {
    tag 'ena_validation_summary'
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path records, stageAs: "records/*"

    output:
    path 'ena_validation_results_mqc.tsv', emit: multiqc
    path 'ena_run_summary.tsv', emit: run_summary
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    collate_ena_validation.py summary \\
        --input 'records/*.ena_validation_result.tsv' \\
        --output ena_validation_results_mqc.tsv \\
        --run-summary ena_run_summary.tsv

    printf '"%s":\n    python: "%s"\n    collate_ena_validation: "1.0.0"\n' \\
        "${task.process}" "\$(python --version | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    printf 'assembly_prefix\tog_id\tvalidation_mode\tena_study\ttable2asn_status\tconversion_status\tpreflight_status\twebin_status\twebin_reason\tsubmission_ready\tvalidation_attempt\n' > ena_validation_results_mqc.tsv
    cp ena_validation_results_mqc.tsv ena_run_summary.tsv
    printf '"%s":\n    python: "stub"\n    collate_ena_validation: "stub"\n' "${task.process}" > versions.yml
    """
}
