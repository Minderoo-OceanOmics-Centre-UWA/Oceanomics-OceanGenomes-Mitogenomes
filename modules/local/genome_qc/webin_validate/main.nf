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

    # A hung ENA request used to sit until Slurm killed the whole task: in
    # batch-02 one OG16 call burned its entire 4 h walltime, and the automatic
    # task retry then passed in 4 s. So bound each call and retry transient
    # failures here, where a second attempt costs seconds rather than hours.
    #
    # The container ships busybox timeout, not GNU coreutils (verified against
    # quay.io/biocontainers/ena-webin-cli:9.0.3): positional seconds only, no
    # --signal/--kill-after long options, and a timeout exits 143 (128+SIGTERM)
    # rather than GNU's 124. Both codes are matched so the conda path, which
    # does supply GNU timeout, behaves identically.
    webin_timeout_s=${params.webin_validate_timeout_seconds}
    max_attempts=${params.webin_validate_max_attempts}
    attempts=0
    webin_rc=0
    webin_status="PASS"
    reason="validated"

    while [ "\${attempts}" -lt "\${max_attempts}" ]; do
        attempts=\$((attempts + 1))

        # Reset per-attempt state: the classifier below greps webin_output, so a
        # report left behind by a failed attempt would misclassify a later one.
        rm -rf webin_output
        mkdir -p webin_output
        : > "\${log_file}"

        timeout -k 30 "\${webin_timeout_s}" ena-webin-cli \\
            -context sequence \\
            -manifest "\${manifest}" \\
            -inputDir . \\
            -outputDir webin_output \\
            -userName "\$WEBIN_USERNAME" \\
            -passwordEnv WEBIN_PASSWORD \\
            -validate > "\${log_file}" 2>&1
        webin_rc=\$?

        if [ "\${webin_rc}" -eq 0 ]; then
            webin_status="PASS"
            reason="validated"
            break
        fi

        if [ "\${webin_rc}" -eq 143 ] || [ "\${webin_rc}" -eq 124 ]; then
            webin_status="FAIL_INFRASTRUCTURE"
            reason="webin_timeout"
        elif grep -RqiE 'validation[^[:alnum:]]*(error|fail)|(^|[^[:alpha:]])(ERROR|INVALID)([^[:alpha:]]|\$)' webin_output "\${log_file}" 2>/dev/null; then
            # Deterministic: the flatfile itself is rejected, so retrying the
            # same input cannot change the answer.
            webin_status="FAIL_WEBIN"
            reason="validation_failed"
            break
        else
            webin_status="FAIL_INFRASTRUCTURE"
            reason="webin_or_network_failure"
        fi

        if [ "\${attempts}" -lt "\${max_attempts}" ]; then
            sleep \$((30 * attempts))
        fi
    done

    # Attempt count goes in the log, not the status TSV: the TSV schema is read
    # by bin/collate_ena_validation.py and is deliberately left unchanged.
    printf 'webin-cli attempts: %s (last exit %s)\n' "\${attempts}" "\${webin_rc}" >> "\${log_file}"

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
