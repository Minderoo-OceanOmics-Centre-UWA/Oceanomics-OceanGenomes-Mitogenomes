process WEBIN_VALIDATE {
    container 'quay.io/biocontainers/ena-webin-cli:9.0.3--hdfd78af_0'

    tag "$meta.id"
    label 'process_low'
    conda 'bioconda::ena-webin-cli=9.0.3'
    secret 'WEBIN_USERNAME'
    secret 'WEBIN_PASSWORD'

    input:
    tuple val(meta), path(embl_file)
    val validation_attempt

    output:
    tuple val(meta), path("validated/*.embl.gz"), optional: true, emit: validated_flatfile
    tuple val(meta), path("*.webin_manifest.txt"), emit: manifest
    tuple val(meta), path("*.webin_status.tsv"), emit: status
    tuple val(meta), path("*.webin_validate.log"), emit: log
    tuple val(meta), path("webin_output"), emit: reports
    tuple val(meta), path("21_webin_validate.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Stem on the full seq id so the manifest, status and log sit alongside
    // the <full_seqid>.embl.gz that ENA_FLATFILE writes.
    def prefix = meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id
    // Per-technology child study, resolved from the candidate by EnaTargets.
    def ena_study = meta.ena_study?.toString()?.trim()
    if (!ena_study) {
        error "meta.ena_study is not set for ${meta.id}: cannot build a Webin manifest"
    }
    """
    set +e
    prefix="${prefix}"
    manifest="\${prefix}.webin_manifest.txt"
    status_file="\${prefix}.webin_status.tsv"
    log_file="\${prefix}.webin_validate.log"
    mkdir -p webin_output validated

    {
        printf 'STUDY\t%s\n' "${ena_study}"
        printf 'NAME\t%s\n' "\${prefix}"
        printf 'DESCRIPTION\tOceanOmics mitochondrial genome %s\n' "\${prefix}"
        printf 'FLATFILE\t%s\n' "${embl_file.name}"
    } > "\${manifest}"

    ena-webin-cli \\
        -context sequence \\
        -manifest "\${manifest}" \\
        -inputDir . \\
        -outputDir webin_output \\
        -userName "\$WEBIN_USERNAME" \\
        -passwordEnv WEBIN_PASSWORD \\
        -validate > "\${log_file}" 2>&1
    webin_rc=\$?

    webin_status="PASS"
    reason="validated"
    if [ "\${webin_rc}" -ne 0 ]; then
        if grep -RqiE 'validation[^[:alnum:]]*(error|fail)|(^|[^[:alpha:]])(ERROR|INVALID)([^[:alpha:]]|\$)' webin_output "\${log_file}" 2>/dev/null; then
            webin_status="FAIL_WEBIN"
            reason="validation_failed"
        else
            webin_status="FAIL_INFRASTRUCTURE"
            reason="webin_or_network_failure"
        fi
    fi

    if [ "\${webin_status}" = "PASS" ]; then
        cp "$embl_file" validated/
    fi
    # The package build copies these counts into package_metadata.json.
    error_count=\$(grep -RihE '(^|[^[:alpha:]])(ERROR|INVALID)([^[:alpha:]]|\$)' webin_output "\${log_file}" 2>/dev/null | wc -l)
    warning_count=\$(grep -RihE '(^|[^[:alpha:]])WARNING([^[:alpha:]]|\$)' webin_output "\${log_file}" 2>/dev/null | wc -l)
    printf 'sample\tstatus\treason\twebin_exit\tvalidation_attempt\terror_count\twarning_count\twebin_cli_version\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t9.0.3\n' \\
        "\${prefix}" "\${webin_status}" "\${reason}" "\${webin_rc}" "${validation_attempt}" \\
        "\${error_count}" "\${warning_count}" > "\${status_file}"

    printf '%s\n' '<tr><td>ENA Webin validation</td><td><samp>ena-webin-cli -context sequence -validate</samp></td><td>Validates the ENA flat file without submitting it for ${meta.id}.</td></tr>' > 21_webin_validate.tool_params_mqcrow.html
    printf '"%s":\n    webin-cli: "9.0.3"\n' "${task.process}" > versions.yml
    exit 0
    """

    stub:
    def prefix = meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id ?: 'stub'
    """
    mkdir -p validated webin_output
    cp "$embl_file" validated/
    printf 'STUDY\t%s\nNAME\t%s\nFLATFILE\t%s\n' "${meta.ena_study}" "${prefix}" "${embl_file.name}" > ${prefix}.webin_manifest.txt
    printf 'sample\tstatus\treason\twebin_exit\tvalidation_attempt\terror_count\twarning_count\twebin_cli_version\n%s\tPASS\tvalidated\t0\t%s\t0\t0\tstub\n' "${prefix}" "${validation_attempt}" > ${prefix}.webin_status.tsv
    printf 'Stub Webin validation passed\n' > ${prefix}.webin_validate.log
    printf 'PASS\n' > webin_output/validation.txt
    printf '%s\n' '<tr><td>ENA Webin validation</td><td><samp>stub</samp></td><td>Stub Webin validation for ${meta.id}.</td></tr>' > 21_webin_validate.tool_params_mqcrow.html
    printf '"%s":\n    webin-cli: "stub"\n' "${task.process}" > versions.yml
    """
}
