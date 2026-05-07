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
    tuple val(meta), path("13_push_lca_blast_results.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                   , emit: versions

    script:
    def effective_args = "${config} ${meta.id} ${lca_results} ${blast_results}"
    """
    # Push the results to SQL database
    push_lca_blast_results.py \\
        $config \\
        ${meta.id} \\
        ${lca_results} \\
        ${blast_results} \\
        > ${meta.id}.lca_blast.upload.txt

    cat <<-END_TOOL_PARAMS > 13_push_lca_blast_results.tool_params_mqcrow.html
    <tr><td>Push LCA BLAST Results</td><td><samp>${effective_args}</samp></td><td>Uploads combined LCA and filtered BLAST results for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control these upload scripts.
    END_VERSIONS
    """
    
    stub:
    def effective_args = "${config} ${meta.id} ${lca_results} ${blast_results}"
    """
    : > ${meta.id}.lca_blast.upload.txt
    cat <<-END_TOOL_PARAMS > 13_push_lca_blast_results.tool_params_mqcrow.html
    <tr><td>Push LCA BLAST Results</td><td><samp>${effective_args}</samp></td><td>Uploads combined LCA and filtered BLAST results for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upload: "stub"
    END_VERSIONS
    """
    }


    
