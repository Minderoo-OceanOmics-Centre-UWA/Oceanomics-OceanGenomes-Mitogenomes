process EXTRACT_GENES_GFF {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(gff)


    output:
    tuple val(meta), path('genes'), emit: genes_dir
    path "versions.yml"                                                , emit: versions

    script:
    def asm = (meta.mt_assembly_prefix ?: meta.sample_id ?: fasta.baseName)

    """
    extract_genes_gff.py \\
        --fasta ${fasta} \\
        --gff ${gff} \\
        --outdir . \\
        --assembly ${asm}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}