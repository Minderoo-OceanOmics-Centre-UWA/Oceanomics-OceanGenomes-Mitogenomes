process ENA_EMBL_PREFLIGHT {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(embl_file)

    output:
    tuple val(meta), path("preflight/*.embl.gz"), optional: true, emit: embl_file
    tuple val(meta), path("*.ena_preflight_status.tsv"), emit: status
    tuple val(meta), path("*.ena_preflight_check.tsv"), emit: checks
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.mt_assembly_prefix ?: meta.id
    """
    set +e
    prefix="${prefix}"
    status_file="\${prefix}.ena_preflight_status.tsv"
    checks_file="\${prefix}.ena_preflight_check.tsv"
    mkdir -p preflight

    gzip -t "$embl_file" >/dev/null 2>&1
    gzip_rc=\$?
    if [ "\${gzip_rc}" -eq 0 ]; then
        id_records=\$(gzip -cd "$embl_file" | grep -c '^ID[[:space:]]' || true)
        sq_records=\$(gzip -cd "$embl_file" | grep -c '^SQ[[:space:]]' || true)
        terminators=\$(gzip -cd "$embl_file" | grep -c '^//\$' || true)
        source_features=\$(gzip -cd "$embl_file" | grep -c '^FT   source[[:space:]]' || true)
        organisms=\$(gzip -cd "$embl_file" | grep -c '/organism=' || true)
        sequence_length=\$(gzip -cd "$embl_file" | awk '/^SQ[[:space:]]/{in_seq=1;next} /^\\/\\//{in_seq=0} in_seq {gsub(/[^A-Za-z]/,""); n+=length(\$0)} END{print n+0}')
    else
        id_records=0
        sq_records=0
        terminators=0
        source_features=0
        organisms=0
        sequence_length=0
    fi

    preflight_status="PASS"
    reason="structure_ok"
    if [ "\${gzip_rc}" -ne 0 ]; then
        preflight_status="FAIL_PREFLIGHT"
        reason="invalid_gzip"
    elif [ "\${id_records}" -lt 1 ] || [ "\${id_records}" -ne "\${sq_records}" ] || [ "\${id_records}" -ne "\${terminators}" ]; then
        preflight_status="FAIL_PREFLIGHT"
        reason="record_structure_mismatch"
    elif [ "\${source_features}" -lt "\${id_records}" ] || [ "\${organisms}" -lt "\${id_records}" ]; then
        preflight_status="FAIL_PREFLIGHT"
        reason="required_annotation_missing"
    elif [ "\${sequence_length}" -le 0 ]; then
        preflight_status="FAIL_PREFLIGHT"
        reason="sequence_missing"
    fi

    printf 'sample\tstatus\treason\n%s\t%s\t%s\n' "\${prefix}" "\${preflight_status}" "\${reason}" > "\${status_file}"
    printf 'sample\tgzip_exit\tid_records\tsq_records\tterminators\tsource_features\torganisms\tsequence_length\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \\
        "\${prefix}" "\${gzip_rc}" "\${id_records}" "\${sq_records}" "\${terminators}" "\${source_features}" "\${organisms}" "\${sequence_length}" > "\${checks_file}"

    if [ "\${preflight_status}" = "PASS" ]; then
        cp "$embl_file" preflight/
    fi
    printf '"%s":\n    gzip: "%s"\n' "${task.process}" "\$(gzip --version | awk 'NR==1{print \$2}')" > versions.yml
    exit 0
    """

    stub:
    def prefix = meta.mt_assembly_prefix ?: meta.id ?: 'stub'
    """
    mkdir -p preflight
    cp "$embl_file" preflight/
    printf 'sample\tstatus\treason\n%s\tPASS\tstructure_ok\n' "${prefix}" > ${prefix}.ena_preflight_status.tsv
    printf 'sample\tgzip_exit\tid_records\tsq_records\tterminators\tsource_features\torganisms\tsequence_length\n%s\t0\t1\t1\t1\t1\t1\t1\n' "${prefix}" > ${prefix}.ena_preflight_check.tsv
    printf '"%s":\n    gzip: "stub"\n' "${task.process}" > versions.yml
    """
}
