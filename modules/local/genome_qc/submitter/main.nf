process SUBMITTER {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    script:
    """
    echo "Submitter placeholder for ${task.process}"
    """

    stub:
    """
    echo "Submitter stub for ${task.process}"
    """
}
