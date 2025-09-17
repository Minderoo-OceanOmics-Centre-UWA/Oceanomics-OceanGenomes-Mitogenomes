process TRANSLATE_GENES {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(genes_dir)


    output:
    tuple val(meta), path('proteins'), emit: proteins_dir
    path "versions.yml"                                                , emit: versions

    script:
    """
    translate_genes.py \\
        --input ${genes_dir} \\
        --outdir . \\
        --table ${params.translation_table ?: 2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p proteins
    : > proteins/${meta.id ?: 'stub'}.faa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
