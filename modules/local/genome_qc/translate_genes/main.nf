process TRANSLATE_GENES {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(genes_dir)


    output:
    tuple val(meta), path('proteins'), emit: proteins_dir
    tuple val(meta), path("18_translate_genes.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                                , emit: versions

    script:
    def effective_args = "--input ${genes_dir} --outdir . --table ${params.translation_table ?: 2}"
    """
    translate_genes.py \\
        --input ${genes_dir} \\
        --outdir . \\
        --table ${params.translation_table ?: 2}

    cat <<-END_TOOL_PARAMS > 18_translate_genes.tool_params_mqcrow.html
    <tr><td>Translate Genes</td><td><samp>${effective_args}</samp></td><td>Translates extracted CDS sequences for ${meta.id} using the configured genetic code.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def effective_args = "--input ${genes_dir} --outdir . --table ${params.translation_table ?: 2}"
    """
    mkdir -p proteins
    : > proteins/${meta.id ?: 'stub'}.faa
    cat <<-END_TOOL_PARAMS > 18_translate_genes.tool_params_mqcrow.html
    <tr><td>Translate Genes</td><td><samp>${effective_args}</samp></td><td>Translates extracted CDS sequences for ${meta.id} using the configured genetic code.</td></tr>
    END_TOOL_PARAMS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
