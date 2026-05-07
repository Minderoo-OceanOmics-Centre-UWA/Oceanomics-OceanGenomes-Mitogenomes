process GETORGANELLE_CONFIG {
    tag "${organelle_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    val(organelle_type)

    output:
    tuple val(organelle_type), path("getorganelle"), emit: db
    path "01_getorganelle_config.tool_params_mqcrow.html", emit: tool_params
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def effective_args = [args, "-a ${organelle_type}", "--config-dir getorganelle"].findAll { it?.trim() }.join(' ')
    """
    get_organelle_config.py \\
        $args \\
        -a ${organelle_type} \\
        --config-dir getorganelle

    cat <<-END_TOOL_PARAMS > 01_getorganelle_config.tool_params_mqcrow.html
    <tr><td>GetOrganelle Config</td><td><samp>${effective_args}</samp></td><td>Downloads the ${organelle_type} GetOrganelle database into the local config directory.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def effective_args = [args, "-a ${organelle_type}", "--config-dir getorganelle"].findAll { it?.trim() }.join(' ')
    """
    mkdir -p getorganelle/{LabelDatabase,SeedDatabase}
    touch getorganelle/{LabelDatabase,SeedDatabase}/${organelle_type}.fasta

    cat <<-END_TOOL_PARAMS > 01_getorganelle_config.tool_params_mqcrow.html
    <tr><td>GetOrganelle Config</td><td><samp>${effective_args}</samp></td><td>Downloads the ${organelle_type} GetOrganelle database into the local config directory.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """
}
