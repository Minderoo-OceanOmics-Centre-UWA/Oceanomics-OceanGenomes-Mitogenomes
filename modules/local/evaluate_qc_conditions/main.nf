process EVALUATE_QC_CONDITIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(blast_table), path(annotation_csv)

    output:
    tuple val(meta), path("species_name.txt"), path("proceed_qc.txt"), emit: evaluation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    evaluate_qc_conditions.py \\
        --blast-table ${blast_table} \\
        --annotation-csv ${annotation_csv} \\
        --output-species species_name.txt \\
        --output-proceed proceed_qc.txt \\
        --output-versions versions.yml
    """

    stub:
    """
    echo "test_species" > species_name.txt
    echo "true" > proceed_qc.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}