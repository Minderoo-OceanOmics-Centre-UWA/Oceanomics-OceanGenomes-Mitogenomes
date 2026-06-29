process GEN_FILES_TABLE2ASN {
    container 'docker://staphb/ncbi-table2asn:latest'

    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::table2asn"

    input:
    tuple val(meta), path(sample_fa), path(sample_tbl), path(sample_cmt), path(sample_src), val(circular)
    path sample_sbt // generic sbt file for the oceanomics lab that you can generate from GenBank.

    output:
    tuple val(meta), path("*.sqn")    , emit: sqn_file
    // A clean validation produces no .val, so it is an optional output; absence
    // means the sample passed table2asn validation with no findings.
    tuple val(meta), path("*.val")    , optional: true, emit: val_file
    tuple val(meta), path("*.gbf")    , emit: gbf_file
    tuple val(meta), path("*.qc_flags.tsv"), emit: qc_flags
    tuple val(meta), path("19_table2asn.tool_params_mqcrow.html"), emit: tool_params

    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def LOGFILE = "${meta.id}.table2asn.log"
    def sample_out = "${meta.mt_assembly_prefix}.sqn"
    // Drive the topology modifier from the resolved per-sample circularity verdict
    // instead of asserting [topology=circular] on every assembly. A genuinely
    // circular molecule is also a complete one, so [completeness=complete] is added
    // to satisfy the SEQ_INST.CompleteCircleProblem validator check; anything not
    // confirmed circular is declared linear (no false circular claim, no warning).
    def topology_mod = (circular?.toString()?.trim() == 'true') ? '[topology=circular] [completeness=complete]' : '[topology=linear]'
    // Per-sample mitochondrial genetic code (meta.genetic_code: 2=vertebrate,
    // 4=coral/coelenterate, 9=echinoderm). Must match the [mgcode] in the FASTA
    // header set by FORMAT_FILES so CDS are translated/validated under the right code.
    def mgcode = meta.genetic_code ?: 2
    def effective_args = "-indir . -euk -J -t ${sample_sbt} -i ${sample_fa} -f ${sample_tbl} -w ${sample_cmt} -src-file ${sample_src} -o ${sample_out} -M n -j '[mgcode=${mgcode}] [location=mitochondrion] ${topology_mod}' -V vb -Z -W"

    """
    # Run table2asn with correct syntax
    table2asn \\
        -indir . \\
        -euk \\
        -J \\
        -t "$sample_sbt" \\
        -i "$sample_fa" \\
        -f "$sample_tbl" \\
        -w "$sample_cmt" \\
        -src-file "$sample_src" \\
        -o "$sample_out" \\
        -M n \\
        -j "[mgcode=${mgcode}] [location=mitochondrion] ${topology_mod}" \\
        -V vb \\
        -Z \\
        -W \\

    # Flag missing-stop-codon (SEQ_FEAT.NoStop) errors from the validation report
    # for manual investigation. These are not blocked here -- they are surfaced in a
    # per-sample TSV (and MultiQC) so the CDS can be reviewed/fixed after the run.
    # A clean validation writes no .val, which is a success, so detect the report
    # with a glob loop that cannot trip 'set -e -o pipefail' when none exists.
    val_file=""
    for f in *.val; do [ -e "\$f" ] && { val_file="\$f"; break; }; done
    flags="${meta.mt_assembly_prefix}.qc_flags.tsv"
    nostop=0
    feats=""
    if [ -n "\${val_file}" ] && [ -s "\${val_file}" ]; then
        nostop=\$(grep -c 'SEQ_FEAT.NoStop' "\${val_file}" || true)
        feats=\$(grep 'SEQ_FEAT.NoStop' "\${val_file}" | sed -E 's/.*CDS: //; s/ *[<[].*//' | paste -sd';' - || true)
    fi
    {
        printf 'sample\\tcircular\\tnostop_count\\tnostop_features\\n'
        printf '%s\\t%s\\t%s\\t%s\\n' "${meta.mt_assembly_prefix}" "${circular}" "\${nostop:-0}" "\${feats}"
    } > "\${flags}"

    cat <<-END_TOOL_PARAMS > 19_table2asn.tool_params_mqcrow.html
    <tr><td>Table2ASN</td><td><samp>${effective_args}</samp></td><td>Generates GenBank submission files and validation reports for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
    
    stub:
    def sample_out = "${meta.mt_assembly_prefix ?: (meta.id ?: 'stub')}.sqn"
    def topology_mod = (circular?.toString()?.trim() == 'true') ? '[topology=circular] [completeness=complete]' : '[topology=linear]'
    def mgcode = meta.genetic_code ?: 2
    def effective_args = "-indir . -euk -J -t ${sample_sbt} -i ${sample_fa} -f ${sample_tbl} -w ${sample_cmt} -src-file ${sample_src} -o ${sample_out} -M n -j '[mgcode=${mgcode}] [location=mitochondrion] ${topology_mod}' -V vb -Z -W"
    """
    prefix=${meta.mt_assembly_prefix ?: (meta.id ?: "stub")}
    : > \${prefix}.sqn
    : > \${prefix}.val
    : > \${prefix}.gbf
    printf 'sample\\tcircular\\tnostop_count\\tnostop_features\\n%s\\t%s\\t0\\t\\n' "\${prefix}" "${circular}" > \${prefix}.qc_flags.tsv
    cat <<-END_TOOL_PARAMS > 19_table2asn.tool_params_mqcrow.html
    <tr><td>Table2ASN</td><td><samp>${effective_args}</samp></td><td>Generates GenBank submission files and validation reports for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        table2asn: "stub"
    END_VERSIONS
    """
}
