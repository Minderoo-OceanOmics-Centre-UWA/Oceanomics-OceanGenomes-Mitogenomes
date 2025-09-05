process PUSH_LCA_BLAST_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(lca_results), path(blast_results)
    path config

    output:
    path "${meta.id}.lca_blast.upload.txt", emit: upload
    path "versions.yml"                   , emit: versions

    script:
    """
    # Push the results to SQL database
    push_lca_blast_results.py \\
        $config \\
        ${meta.id} \\
        ${lca_results} \\
        ${blast_results} \\
        > ${meta.id}.lca_blast.upload.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control these upload scripts.
    END_VERSIONS
    """
    }


    