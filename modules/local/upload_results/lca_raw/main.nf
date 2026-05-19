process PUSH_LCA_RAW_RESULTS {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9 conda-forge::psycopg2 conda-forge::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(lca_raw_files)
    path config

    output:
    path "${meta.id}.lca_raw.upload.txt", emit: upload
    tuple val(meta), path("15_push_lca_raw_results.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def files_str = lca_raw_files instanceof List ? lca_raw_files.collect{ it.toString() }.join(' ') : lca_raw_files.toString()
    def effective_args = [args, config, meta.id, files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    push_lca_raw_results.py \\
        $args \\
        $config \\
        ${meta.id} \\
        ${files_str} \\
        > ${meta.id}.lca_raw.upload.txt

    cat <<-END_TOOL_PARAMS > 15_push_lca_raw_results.tool_params_mqcrow.html
    <tr><td>Push LCA raw results</td><td><samp>${effective_args}</samp></td><td>Uploads raw per-hit LCA rows for ${meta.id} to the lca_raw_results table.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        push_lca_raw_results: "1.0.0"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def files_str = lca_raw_files instanceof List ? lca_raw_files.collect{ it.toString() }.join(' ') : lca_raw_files.toString()
    def effective_args = [args, config, meta.id, files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    touch ${meta.id}.lca_raw.upload.txt

    cat <<-END_TOOL_PARAMS > 15_push_lca_raw_results.tool_params_mqcrow.html
    <tr><td>Push LCA raw results</td><td><samp>${effective_args}</samp></td><td>Uploads raw per-hit LCA rows for ${meta.id} to the lca_raw_results table.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        push_lca_raw_results: "1.0.0"
    END_VERSIONS
    """
}
