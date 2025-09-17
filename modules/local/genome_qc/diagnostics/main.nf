process DIAGNOSTICS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    script:
    """
    echo "Diagnostics placeholder for ${task.process}"
    """

    stub:
    """
    echo "Diagnostics stub for ${task.process}"
    """
}
