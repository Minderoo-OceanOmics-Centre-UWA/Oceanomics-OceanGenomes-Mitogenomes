process PREPARE_LCA_DATABASES {
    tag 'taxonomy-cache'
    label 'process_single'
    label 'error_retry'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0' :
        'quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0' }"

    // This process validates the persistent shared cache on every workflow run.
    // The lock makes concurrent runs safe; existing valid files are not downloaded.
    cache false

    input:
    val cache_dir
    path prepare_script

    output:
    path 'lca_cache_manifest.json', emit: manifest
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${prepare_script} \\
        --cache-dir '${cache_dir}' \\
        --manifest lca_cache_manifest.json

    printf '"%s":\n    python: "%s"\n' \
        '${task.process}' \
        "\$(python --version | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf '{"cache_dir":"%s","cache_schema_version":1,"cache_signature":"0000000000000000000000000000000000000000000000000000000000000000","files":{},"status":"stub"}\n' '${cache_dir}' > lca_cache_manifest.json
    printf '"%s":\n    python: "stub"\n' '${task.process}' > versions.yml
    """
}
