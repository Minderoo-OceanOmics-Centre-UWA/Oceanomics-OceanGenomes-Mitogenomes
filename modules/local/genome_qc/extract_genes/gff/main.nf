process EXTRACT_GENES_GFF {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(gff)


    output:
    tuple val(meta), path('genes'), emit: genes_dir
    tuple val(meta), path("17_extract_genes_gff.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                                , emit: versions

    script:
    def asm = (meta.mt_assembly_prefix ?: meta.sample_id ?: fasta.baseName)
    def effective_args = "--fasta ${fasta} --gff ${gff} --outdir . --assembly ${asm}"

    """
    extract_genes_gff.py \\
        --fasta ${fasta} \\
        --gff ${gff} \\
        --outdir . \\
        --assembly ${asm}

    cat <<-END_TOOL_PARAMS > 17_extract_genes_gff.tool_params_mqcrow.html
    <tr><td>Extract Genes GFF</td><td><samp>${effective_args}</samp></td><td>Extracts gene sequences from GFF annotation for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def asm = (meta.mt_assembly_prefix ?: meta.sample_id ?: fasta.baseName)
    def effective_args = "--fasta ${fasta} --gff ${gff} --outdir . --assembly ${asm}"
    """
    mkdir -p genes
    cat <<-END_TOOL_PARAMS > 17_extract_genes_gff.tool_params_mqcrow.html
    <tr><td>Extract Genes GFF</td><td><samp>${effective_args}</samp></td><td>Extracts gene sequences from GFF annotation for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    : > versions.yml
    """
}
