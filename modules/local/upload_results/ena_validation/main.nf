process PUSH_ENA_VALIDATION_RESULTS {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(validation_record)
    path config

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.ena_validation.upload.txt"), emit: upload
    tuple val(meta), path("15_push_ena_validation_results.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set +e
    push_ena_validation_results.py \\
        "$config" \\
        "$validation_record" \\
        > '${meta.mt_assembly_prefix}.ena_validation.upload.txt' 2>&1
    upload_rc=\$?
    printf 'UPLOAD_EXIT=%s\n' "\${upload_rc}" >> '${meta.mt_assembly_prefix}.ena_validation.upload.txt'

    printf '%s\n' '<tr><td>Push ENA Validation Results</td><td><samp>push_ena_validation_results.py &lt;config&gt; &lt;record&gt;</samp></td><td>Uploads normalized ENA validation history for ${meta.id}.</td></tr>' > 15_push_ena_validation_results.tool_params_mqcrow.html
    printf '"%s":\n    python: "%s"\n    push_ena_validation_results: "1.0.0"\n' \\
        "${task.process}" "\$(python --version | awk '{print \$2}')" > versions.yml
    exit 0
    """

    stub:
    """
    printf '✅ Success: inserted ENA validation attempt for %s\nUPLOAD_EXIT=0\n' '${meta.mt_assembly_prefix}' > '${meta.mt_assembly_prefix}.ena_validation.upload.txt'
    printf '%s\n' '<tr><td>Push ENA Validation Results</td><td><samp>stub</samp></td><td>Stub ENA validation upload for ${meta.id}.</td></tr>' > 15_push_ena_validation_results.tool_params_mqcrow.html
    printf '"%s":\n    python: "stub"\n    push_ena_validation_results: "stub"\n' "${task.process}" > versions.yml
    """
}
