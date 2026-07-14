process ENA_FLATFILE {
    container 'quay.io/biocontainers/emboss:6.6.0--hd9b00e3_12'

    tag "$meta.id"
    label 'process_low'
    conda 'bioconda::emboss=6.6.0'

    input:
    tuple val(meta), path(gbf)

    output:
    tuple val(meta), path("*.embl.gz"), optional: true, emit: embl_file
    tuple val(meta), path("*.ena_conversion_status.tsv"), emit: status
    tuple val(meta), path("*.ena_conversion_check.tsv"), emit: checks
    tuple val(meta), path("*.ena_conversion.log"), emit: log
    tuple val(meta), path("20_ena_flatfile.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.mt_assembly_prefix ?: meta.id
    """
    set +e
    prefix="${prefix}"
    raw="\${prefix}.embl"
    compressed="\${raw}.gz"
    log="\${prefix}.ena_conversion.log"
    status_file="\${prefix}.ena_conversion_status.tsv"
    checks_file="\${prefix}.ena_conversion_check.tsv"

    seqret -auto -feature -sformat genbank -osformat embl -sequence "$gbf" -outseq "\${raw}" > "\${log}" 2>&1
    seqret_rc=\$?

    input_records=\$(grep -c '^LOCUS[[:space:]]' "$gbf" || true)
    output_records=\$(grep -c '^ID[[:space:]]' "\${raw}" 2>/dev/null || true)
    terminators=\$(grep -c '^//\$' "\${raw}" 2>/dev/null || true)
    input_length=\$(awk '/^ORIGIN/{in_seq=1;next} /^\\/\\//{in_seq=0} in_seq {gsub(/[^A-Za-z]/,""); n+=length(\$0)} END{print n+0}' "$gbf")
    output_length=\$(awk '/^SQ[[:space:]]/{in_seq=1;next} /^\\/\\//{in_seq=0} in_seq {gsub(/[^A-Za-z]/,""); n+=length(\$0)} END{print n+0}' "\${raw}" 2>/dev/null)
    input_features=\$(grep -Ec '^     [A-Za-z_][A-Za-z_0-9]*[[:space:]]' "$gbf" || true)
    output_features=\$(grep -Ec '^FT   [^[:space:]]+[[:space:]]' "\${raw}" 2>/dev/null || true)

    missing_qualifiers=""
    for qualifier in organism mol_type organelle gene product transl_table codon_start translation; do
        if grep -q "/\${qualifier}=" "$gbf" && ! grep -q "/\${qualifier}=" "\${raw}"; then
            missing_qualifiers="\${missing_qualifiers}\${missing_qualifiers:+,}\${qualifier}"
        fi
    done

    source_ok=0
    organism_ok=0
    topology_ok=1
    grep -Eq '^FT   source[[:space:]]' "\${raw}" 2>/dev/null && source_ok=1
    grep -q '/organism=' "\${raw}" 2>/dev/null && organism_ok=1
    if grep -Eq '^LOCUS.*[[:space:]]circular[[:space:]]' "$gbf" && ! grep -Eq '^ID.*; circular;' "\${raw}" 2>/dev/null; then
        topology_ok=0
    fi

    conversion_status="PASS"
    reason="ok"
    if [ "\${seqret_rc}" -ne 0 ]; then
        conversion_status="FAIL_CONVERSION"
        reason="seqret_exit_\${seqret_rc}"
    elif [ "\${input_records}" -lt 1 ] || [ "\${input_records}" -ne "\${output_records}" ] || [ "\${output_records}" -ne "\${terminators}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="record_structure_mismatch"
    elif [ "\${input_length}" -le 0 ] || [ "\${input_length}" -ne "\${output_length:-0}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="sequence_length_mismatch"
    elif [ "\${input_features}" -ne "\${output_features}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="feature_count_mismatch"
    elif [ "\${source_ok}" -ne 1 ] || [ "\${organism_ok}" -ne 1 ] || [ "\${topology_ok}" -ne 1 ] || [ -n "\${missing_qualifiers}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="required_annotation_not_preserved"
    fi

    {
        printf 'sample\tstatus\treason\tseqret_exit\n'
        printf '%s\t%s\t%s\t%s\n' "\${prefix}" "\${conversion_status}" "\${reason}" "\${seqret_rc}"
    } > "\${status_file}"
    {
        printf 'sample\tinput_records\toutput_records\tterminators\tinput_length\toutput_length\tinput_features\toutput_features\tsource_ok\torganism_ok\ttopology_ok\tmissing_qualifiers\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "\${prefix}" "\${input_records}" "\${output_records}" "\${terminators}" "\${input_length}" "\${output_length:-0}" "\${input_features}" "\${output_features}" "\${source_ok}" "\${organism_ok}" "\${topology_ok}" "\${missing_qualifiers}"
    } > "\${checks_file}"

    if [ "\${conversion_status}" = "PASS" ]; then
        gzip -n "\${raw}"
    else
        rm -f "\${raw}" "\${compressed}"
    fi

    printf '%s\n' '<tr><td>ENA flat-file conversion</td><td><samp>seqret -feature -sformat genbank -osformat embl</samp></td><td>Converts the validated table2asn GenBank flat file for ENA and checks annotation preservation for ${meta.id}.</td></tr>' > 20_ena_flatfile.tool_params_mqcrow.html
    emboss_version=\$(seqret -version 2>&1 | awk 'NR==1{print \$NF}')
    printf '"%s":\n    emboss: "%s"\n' "${task.process}" "\${emboss_version}" > versions.yml
    exit 0
    """

    stub:
    def prefix = meta.mt_assembly_prefix ?: meta.id ?: 'stub'
    """
    : | gzip -n > ${prefix}.embl.gz
    printf 'sample\tstatus\treason\tseqret_exit\n%s\tPASS\tok\t0\n' "${prefix}" > ${prefix}.ena_conversion_status.tsv
    printf 'sample\tinput_records\toutput_records\tterminators\tinput_length\toutput_length\tinput_features\toutput_features\tsource_ok\torganism_ok\ttopology_ok\tmissing_qualifiers\n%s\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t\n' "${prefix}" > ${prefix}.ena_conversion_check.tsv
    printf 'Stub EMBOSS conversion passed\n' > ${prefix}.ena_conversion.log
    printf '%s\n' '<tr><td>ENA flat-file conversion</td><td><samp>stub</samp></td><td>Stub ENA conversion for ${meta.id}.</td></tr>' > 20_ena_flatfile.tool_params_mqcrow.html
    printf '"%s":\n    emboss: "stub"\n' "${task.process}" > versions.yml
    """
}
