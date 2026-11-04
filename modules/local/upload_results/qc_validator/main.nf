process PUSH_QC_VALIDATOR {
    tag "$meta.id"
    label 'process_upload'

    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(validation_record)
    path config

    output:
    tuple val(meta), path("${meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id}.qc_validator.upload.txt"), emit: upload
    tuple val(meta), path("16_push_qc_validator.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def full_seqid = meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id
    """
    set +e
    push_qc_validator.py \\
        "$config" \\
        "$validation_record" \\
        > '${full_seqid}.qc_validator.upload.txt' 2>&1
    upload_rc=\$?
    printf 'UPLOAD_EXIT=%s\n' "\${upload_rc}" >> '${full_seqid}.qc_validator.upload.txt'

    printf '%s\n' '<tr><td>Push QC Validator</td><td><samp>push_qc_validator.py &lt;config&gt; &lt;record&gt;</samp></td><td>Records the pipeline as lca_validation.validator_2 for ${meta.id} when it cleared every QC gate.</td></tr>' > 16_push_qc_validator.tool_params_mqcrow.html
    printf '"%s":\n    python: "%s"\n    push_qc_validator: "1.0.0"\n' \\
        "${task.process}" "\$(python --version | awk '{print \$2}')" > versions.yml
    # Exit with the script's own code. A failed upload must fail the task: the upload
    # label's errorStrategy retries it once, then ignores it and records it in
    # FAILED_UPLOADS.tsv, so the run carries on and -resume re-runs it. The log above is
    # written either way.
    exit \${upload_rc}
    """

    stub:
    def full_seqid = meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id
    """
    printf "✅ Success: lca_validation validator_2 set to 'QCd-nf-core' for %s\nUPLOAD_EXIT=0\n" '${full_seqid}' > '${full_seqid}.qc_validator.upload.txt'
    printf '%s\n' '<tr><td>Push QC Validator</td><td><samp>stub</samp></td><td>Stub validator_2 upload for ${meta.id}.</td></tr>' > 16_push_qc_validator.tool_params_mqcrow.html
    printf '"%s":\n    python: "stub"\n    push_qc_validator: "stub"\n' "${task.process}" > versions.yml
    """
}
