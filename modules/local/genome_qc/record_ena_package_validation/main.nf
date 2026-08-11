process RECORD_ENA_PACKAGE_VALIDATION {
    tag "$meta.full_seqid:$service"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), val(service), path(status)
    path db_config

    output:
    tuple val(meta), val(service), path("*.ena_validation_recorded.tsv"), emit: receipt
    path "versions.yml", emit: versions

    script:
    """
    record_ena_package_validation.py \
        --config '${db_config}' \
        --status '${status}' \
        --output '${meta.full_seqid}.${service}.ena_validation_recorded.tsv'
    printf '"%s":\n    python: "%s"\n    record_ena_package_validation: "1.0.0"\n' \
        "${task.process}" "\$(python --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf 'full_seqid\tservice\tstatus\trecorded\n%s\t%s\tPASS\ttrue\n' \
        '${meta.full_seqid}' '${service}' > '${meta.full_seqid}.${service}.ena_validation_recorded.tsv'
    printf '"%s":\n    python: "stub"\n    record_ena_package_validation: "stub"\n' "${task.process}" > versions.yml
    """
}
