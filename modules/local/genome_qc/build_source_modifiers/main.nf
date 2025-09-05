process BUILD_SOURCE_MODIFIERS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    val meta
    path config

    output:
    tuple val(meta), path("${meta.id}.bankit_metadata.csv")                , emit: bankit_metadata
    tuple val(meta), path("${meta.id}.bankit_metadata_latlon_cleaned.csv") , emit: cleaned_metadata
    tuple val(meta), path("${meta.mt_assembly_prefix}*.src")                , emit: src_file, optional: true
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    build_source_modifiers.py \\
        --config $config \\
        --og-id "${meta.id}" \\
        --seq-tech "${meta.sequencing_type}" \\
        --assembly-id "${meta.mt_assembly_prefix}" \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output src_files
    : > output/bankit_metadata.csv
    : > output/bankit_metadata_latlon_cleaned.csv
    : > src_files/dummy.src

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //' 2>/dev/null || echo "3.9.0")
    END_VERSIONS
    """
}
