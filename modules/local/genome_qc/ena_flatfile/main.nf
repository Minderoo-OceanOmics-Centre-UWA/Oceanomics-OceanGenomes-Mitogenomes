process ENA_FLATFILE {
    container 'quay.io/biocontainers/emboss:6.6.0--hd9b00e3_12'

    tag "$meta.id"
    label 'process_low'
    conda 'bioconda::emboss=6.6.0'

    input:
    tuple val(meta), path(gbf)

    output:
    tuple val(meta), path("*.embl.gz"), optional: true, emit: embl_file
    tuple val(meta), path("*.embl"), optional: true, emit: embl_file_plain
    tuple val(meta), path("*.ena_conversion_status.tsv"), emit: status
    tuple val(meta), path("*.ena_conversion_check.tsv"), emit: checks
    tuple val(meta), path("*.ena_conversion.log"), emit: log
    tuple val(meta), path("20_ena_flatfile.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id
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
    postprocess_rc=0

    # seqret's GenBank->EMBL conversion has two known gaps for pre-accession
    # organelle submissions, both of which ENA's flatfile validator rejects:
    if [ -e "\${raw}" ]; then
        # 1. It never writes an AC line when the source GenBank ACCESSION field
        # is blank (true for anything not yet accessioned), but ENA requires the
        # AC block to appear exactly once regardless. Insert the standard
        # not-yet-accessioned placeholder right after the ID block.
        # Both lines are required and they are not interchangeable: "AC   ;" is
        # the empty AC block webin counts, while "AC * _<entry>" is the separate
        # entry-name directive. Emitting only the directive leaves the record
        # with zero AC blocks and webin fails with
        # "ERROR: Block AC must occur exactly once".
        awk -v entry="\${prefix}" '
            /^ID   / {
                sub(/^ID   [^;]+;/, "ID   " entry ";")
                print
                next
            }
            /^AC / {
                if (!seen_ac++) {
                    print "AC   ;"
                    print "XX"
                    print "AC * _" entry
                }
                next
            }
            index(\$0, "/note=\\"*geo_loc_name: ") {
                value = \$0
                sub(/^FT +/, "", value)
                sub("^/note=\\"[*]geo_loc_name: ", "", value)
                while (value !~ /"\$/ && (getline continuation) > 0) {
                    sub(/^FT +/, "", continuation)
                    value = value " " continuation
                }
                sub(/"\$/, "", value)
                print "FT                   /geo_loc_name=\\"" value "\\""
                next
            }
            {print}
            END {
                if (!seen_ac) exit 42
            }
        ' "\${raw}" > "\${raw}.tmp"
        awk_rc=\$?
        if [ "\${awk_rc}" -eq 42 ]; then
            # Sentinel: seqret emitted no AC line at all (the usual case for a
            # blank ACCESSION), so the pass above had nothing to rewrite. Its
            # output is otherwise complete, so insert the placeholder into that
            # result -- re-running against "\${raw}" would silently discard the
            # ID-line and geo_loc_name fixes it just made.
            awk -v entry="\${prefix}" '/^ID   /{print; print "XX"; print "AC   ;"; print "XX"; print "AC * _" entry; next} {print}' "\${raw}.tmp" > "\${raw}.tmp2"
            awk_rc=\$?
            if [ "\${awk_rc}" -eq 0 ]; then
                mv "\${raw}.tmp2" "\${raw}.tmp"
            else
                rm -f "\${raw}.tmp2"
            fi
        fi
        if [ "\${awk_rc}" -eq 0 ]; then
            mv "\${raw}.tmp" "\${raw}"
        else
            # Any other non-zero exit means awk died partway (or before reading
            # anything), leaving a truncated or empty temp file. Keep seqret's
            # output intact and record the real cause, otherwise the structural
            # checks below compare against an empty file and misreport this as a
            # record_structure_mismatch.
            rm -f "\${raw}.tmp"
            postprocess_rc=\${awk_rc}
        fi

        # seqret propagates "SOURCE mitochondrion <organism>" into OS and turns
        # /geo_loc_name into a wrapped "*geo_loc_name" note. The awk pass above
        # restores the qualifier; derive OS/DE from /organism so names containing
        # punctuation never become shell code.
        organism=\$(grep -m1 '/organism=' "\${raw}" | sed -E 's#.*="([^"]+)".*#\\1#')
        awk -v organism="\${organism}" '
            /^OS   / { print "OS   " organism; next }
            /^DE   / { print "DE   " organism " mitochondrion, complete genome"; next }
            { print }
        ' "\${raw}" > "\${raw}.tmp" && mv "\${raw}.tmp" "\${raw}"

        # 2. It defaults the ID line's molecule-type token to "unassigned DNA"
        # instead of deriving it from the source feature, leaving it
        # inconsistent with the preserved /mol_type qualifier. Sync the ID
        # line to whatever /mol_type the source feature actually carries.
        mol_type=\$(grep -m1 '/mol_type=' "\${raw}" | sed -E 's#.*/mol_type="([^"]*)".*#\\1#')
        if [ -n "\${mol_type}" ]; then
            sed -i -E "s/(^ID   [^;]+; SV [0-9]+; [a-z]+;) [^;]+;/\\1 \${mol_type};/" "\${raw}"
        fi
    fi

    input_records=\$(grep -c '^LOCUS[[:space:]]' "$gbf" || true)
    output_records=\$(grep -c '^ID[[:space:]]' "\${raw}" 2>/dev/null || true)
    terminators=\$(grep -c '^//\$' "\${raw}" 2>/dev/null || true)
    input_length=\$(awk '/^ORIGIN/{in_seq=1;next} /^\\/\\//{in_seq=0} in_seq {gsub(/[^A-Za-z]/,""); n+=length(\$0)} END{print n+0}' "$gbf")
    output_length=\$(awk '/^SQ[[:space:]]/{in_seq=1;next} /^\\/\\//{in_seq=0} in_seq {gsub(/[^A-Za-z]/,""); n+=length(\$0)} END{print n+0}' "\${raw}" 2>/dev/null)
    input_features=\$(grep -Ec '^     [A-Za-z_][A-Za-z_0-9]*[[:space:]]' "$gbf" || true)
    output_features=\$(grep -Ec '^FT   [^[:space:]]+[[:space:]]' "\${raw}" 2>/dev/null || true)

    # Compare qualifier COUNTS between the GenBank input and the EMBL output, not
    # mere presence. A presence test passes as long as one feature anywhere in the
    # record still carries the qualifier, so a single CDS losing its /translation
    # is invisible to it -- which is exactly how an ATP6 whose translation began
    # with a gap symbol reached ENA with conversion_status=PASS.
    missing_qualifiers=""
    for qualifier in organism mol_type organelle gene product transl_table codon_start translation geo_loc_name; do
        gbf_count=\$(grep -c "/\${qualifier}=" "$gbf" || true)
        embl_count=\$(grep -c "/\${qualifier}=" "\${raw}" 2>/dev/null || true)
        if [ "\${gbf_count:-0}" -ne "\${embl_count:-0}" ]; then
            missing_qualifiers="\${missing_qualifiers}\${missing_qualifiers:+,}\${qualifier}(\${gbf_count:-0}->\${embl_count:-0})"
        fi
    done

    # seqret never fails on a qualifier it dislikes; it demotes it to free text as
    # /note="*<name>: ...". The post-processing above deliberately restores
    # geo_loc_name, so anything still demoted at this point is a qualifier that
    # would reach ENA as prose and be silently lost (a dropped /translation, an
    # invalid /lat_lon). One catch-all check covers every such qualifier,
    # including ones not in the list above.
    demoted_qualifiers=\$(grep -o '/note="[*][a-z_]*' "\${raw}" 2>/dev/null | sed 's#.*/note="[*]##' | sort -u | paste -sd, -)

    source_ok=0
    organism_ok=0
    topology_ok=1
    grep -Eq '^FT   source[[:space:]]' "\${raw}" 2>/dev/null && source_ok=1
    grep -q '/organism=' "\${raw}" 2>/dev/null && organism_ok=1
    if grep -Eq '^LOCUS.*[[:space:]]circular([[:space:]]|\$)' "$gbf" && ! grep -Eq '^ID.*; circular;' "\${raw}" 2>/dev/null; then
        topology_ok=0
    fi

    conversion_status="PASS"
    reason="ok"
    if [ "\${seqret_rc}" -ne 0 ]; then
        conversion_status="FAIL_CONVERSION"
        reason="seqret_exit_\${seqret_rc}"
    elif [ "\${postprocess_rc}" -ne 0 ]; then
        conversion_status="FAIL_CONVERSION"
        reason="embl_postprocess_exit_\${postprocess_rc}"
    elif [ "\${input_records}" -lt 1 ] || [ "\${input_records}" -ne "\${output_records}" ] || [ "\${output_records}" -ne "\${terminators}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="record_structure_mismatch"
    elif [ "\${input_length}" -le 0 ] || [ "\${input_length}" -ne "\${output_length:-0}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="sequence_length_mismatch"
    elif [ "\${input_features}" -ne "\${output_features}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="feature_count_mismatch"
    elif [ -n "\${demoted_qualifiers}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="qualifier_demoted_by_seqret:\${demoted_qualifiers}"
    elif [ "\${source_ok}" -ne 1 ] || [ "\${organism_ok}" -ne 1 ] || [ "\${topology_ok}" -ne 1 ] || [ -n "\${missing_qualifiers}" ]; then
        conversion_status="FAIL_CONVERSION"
        reason="required_annotation_not_preserved"
    fi

    {
        printf 'sample\tstatus\treason\tseqret_exit\n'
        printf '%s\t%s\t%s\t%s\n' "\${prefix}" "\${conversion_status}" "\${reason}" "\${seqret_rc}"
    } > "\${status_file}"
    {
        printf 'sample\tinput_records\toutput_records\tterminators\tinput_length\toutput_length\tinput_features\toutput_features\tsource_ok\torganism_ok\ttopology_ok\tmissing_qualifiers\tdemoted_qualifiers\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "\${prefix}" "\${input_records}" "\${output_records}" "\${terminators}" "\${input_length}" "\${output_length:-0}" "\${input_features}" "\${output_features}" "\${source_ok}" "\${organism_ok}" "\${topology_ok}" "\${missing_qualifiers}" "\${demoted_qualifiers}"
    } > "\${checks_file}"

    if [ "\${conversion_status}" = "PASS" ]; then
        gzip -n -k "\${raw}"
    else
        rm -f "\${raw}" "\${compressed}"
    fi

    printf '%s\n' '<tr><td>ENA flat-file conversion</td><td><samp>seqret -feature -sformat genbank -osformat embl</samp></td><td>Converts the validated table2asn GenBank flat file for ENA and checks annotation preservation for ${meta.id}.</td></tr>' > 20_ena_flatfile.tool_params_mqcrow.html
    emboss_version=\$(seqret -version 2>&1 | awk 'NR==1{print \$NF}')
    printf '"%s":\n    emboss: "%s"\n' "${task.process}" "\${emboss_version}" > versions.yml
    exit 0
    """

    stub:
    def prefix = meta.full_seqid ?: meta.mt_assembly_prefix ?: meta.id ?: 'stub'
    """
    : > ${prefix}.embl
    gzip -n -k ${prefix}.embl
    printf 'sample\tstatus\treason\tseqret_exit\n%s\tPASS\tok\t0\n' "${prefix}" > ${prefix}.ena_conversion_status.tsv
    printf 'sample\tinput_records\toutput_records\tterminators\tinput_length\toutput_length\tinput_features\toutput_features\tsource_ok\torganism_ok\ttopology_ok\tmissing_qualifiers\tdemoted_qualifiers\n%s\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t\t\n' "${prefix}" > ${prefix}.ena_conversion_check.tsv
    printf 'Stub EMBOSS conversion passed\n' > ${prefix}.ena_conversion.log
    printf '%s\n' '<tr><td>ENA flat-file conversion</td><td><samp>stub</samp></td><td>Stub ENA conversion for ${meta.id}.</td></tr>' > 20_ena_flatfile.tool_params_mqcrow.html
    printf '"%s":\n    emboss: "stub"\n' "${task.process}" > versions.yml
    """
}
