process PUSH_MTDNA_ASSM_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(out_log)
    path config

    output:
    path "${meta.id}.mtdna.upload.txt", emit: upload
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    """
    push_mtdna_assm_results.py \\
        $args \\
        $config \\
        ${mt_assembly_prefix} \\
        $out_log \\
        $fasta \\
        > ${meta.id}.mtdna.upload.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        push_mtdna_assm_results: "1.0.0"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mtdna.upload.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        push_mtdna_assm_results: "1.0.0"
    END_VERSIONS
    """
}