process WEBIN_VALIDATE_GENOME {
    container 'quay.io/biocontainers/ena-webin-cli:9.0.3--hdfd78af_0'

    tag "$meta.full_seqid:$service"
    label 'process_low'
    conda 'bioconda::ena-webin-cli=9.0.3'
    secret 'WEBIN_USERNAME'
    secret 'WEBIN_PASSWORD'

    input:
    tuple val(meta), path(package_dir)
    val service
    val validation_attempt

    output:
    tuple val(meta), val(service), path("*.webin_${service}_status.tsv"), emit: status
    tuple val(meta), val(service), path("webin_${service}_output"), emit: reports
    tuple val(meta), val(service), path("*.webin_${service}.log"), emit: log
    tuple val(meta), path("23_webin_genome_${service}.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def test_arg = service == 'test' ? '-test' : ''
    if (!(service in ['test', 'production'])) {
        error "WEBIN_VALIDATE_GENOME service must be test or production, got '${service}'"
    }
    """
    set +e
    mkdir -p input webin_${service}_output
    cp -L package/* input/
    manifest=\$(find input -maxdepth 1 -name '*.manifest.txt' -type f | head -1)
    status='${meta.full_seqid}.webin_${service}_status.tsv'
    log='${meta.full_seqid}.webin_${service}.log'
    if [ -z "\${manifest}" ]; then
        printf 'full_seqid\tservice\tstatus\treason\twebin_exit\tvalidation_attempt\terror_count\twarning_count\terror_codes\twebin_cli_version\treport_path\n%s\t%s\tFAIL_PACKAGE\tmissing_manifest\t\t%s\t\t\t\t9.0.3\twebin_%s_output\n' \
            '${meta.full_seqid}' '${service}' '${validation_attempt}' '${service}' > "\${status}"
        printf 'Missing package manifest\n' > "\${log}"
        exit 0
    fi

    ena-webin-cli \
        -context genome \
        -manifest "\${manifest}" \
        -inputDir input \
        -outputDir webin_${service}_output \
        -userName "\$WEBIN_USERNAME" \
        -passwordEnv WEBIN_PASSWORD \
        -validate ${test_arg} > "\${log}" 2>&1
    rc=\$?
    result=PASS
    reason=validated
    if [ "\${rc}" -ne 0 ]; then
        if grep -RqiE 'validation[^[:alnum:]]*(error|fail)|(^|[^[:alpha:]])(ERROR|INVALID)([^[:alpha:]]|\$)' webin_${service}_output "\${log}" 2>/dev/null; then
            result=FAIL_WEBIN
            reason=validation_failed
        else
            result=FAIL_INFRASTRUCTURE
            reason=webin_or_network_failure
        fi
    fi
    error_count=\$(grep -RihE '(^|[^[:alpha:]])(ERROR|INVALID)([^[:alpha:]]|\$)' webin_${service}_output "\${log}" 2>/dev/null | wc -l)
    warning_count=\$(grep -RihE '(^|[^[:alpha:]])WARNING([^[:alpha:]]|\$)' webin_${service}_output "\${log}" 2>/dev/null | wc -l)
    printf 'full_seqid\tservice\tstatus\treason\twebin_exit\tvalidation_attempt\terror_count\twarning_count\terror_codes\twebin_cli_version\treport_path\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t9.0.3\t%s\n' \
        '${meta.full_seqid}' '${service}' "\${result}" "\${reason}" "\${rc}" '${validation_attempt}' \
        "\${error_count}" "\${warning_count}" 'webin_${service}_output' > "\${status}"
    printf '%s\n' '<tr><td>ENA Webin ${service}</td><td><samp>ena-webin-cli -context genome -validate ${test_arg}</samp></td><td>Validates the finished genome-context package for ${meta.full_seqid} without submitting it.</td></tr>' > 23_webin_genome_${service}.tool_params_mqcrow.html
    printf '"%s":\n    webin-cli: "9.0.3"\n' "${task.process}" > versions.yml
    exit 0
    """

    stub:
    """
    mkdir -p webin_${service}_output
    printf 'full_seqid\tservice\tstatus\treason\twebin_exit\tvalidation_attempt\terror_count\twarning_count\terror_codes\twebin_cli_version\treport_path\n%s\t%s\tPASS\tvalidated\t0\t%s\t0\t0\t\tstub\t%s\n' \
        '${meta.full_seqid}' '${service}' '${validation_attempt}' 'webin_${service}_output' > '${meta.full_seqid}.webin_${service}_status.tsv'
    printf 'PASS\n' > webin_${service}_output/validation.txt
    printf 'Stub Webin validation\n' > '${meta.full_seqid}.webin_${service}.log'
    printf '%s\n' '<tr><td>ENA Webin ${service}</td><td><samp>stub</samp></td><td>Stub genome-context validation for ${meta.full_seqid}.</td></tr>' > 23_webin_genome_${service}.tool_params_mqcrow.html
    printf '"%s":\n    webin-cli: "stub"\n' "${task.process}" > versions.yml
    """
}
