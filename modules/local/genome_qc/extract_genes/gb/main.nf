process EXTRACT_GENES_GB {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(tbl)


    output:
    tuple val(meta), path('cds'), emit: cds_dir
    tuple val(meta), path('cds/*.fa'), emit: cds_fastas
    tuple val(meta), path('genes'), optional: true, emit: genes_dir
    tuple val(meta), path('genes/*.fa'), optional: true, emit: genes_fastas
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when ?: true


    script:
    def asm = (meta.mt_assembly_prefix ?: meta.sample_id ?: fasta.baseName)

    """
    extract_cds_from_tbl.py \\
        --fasta ${fasta} \\
        --tbl ${tbl} \\
        --outdir . \\
        --assembly ${asm} 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
