process GEN_FILES_TABLE2ASN {
    container 'docker://staphb/ncbi-table2asn:latest'

    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::table2asn"

    input:
    tuple val(meta), path(sample_fa), path(sample_tbl), path(sample_cmt), path(sample_src), val(circular)
    path sample_sbt

    output:
    tuple val(meta), path("*.sqn")    , emit: sqn_file
    tuple val(meta), path("*.val"), emit: val_file
    tuple val(meta), path("*.stats")  , optional: true, emit: stats_file
    tuple val(meta), path("*.dr"), emit: discrepancy_file
    tuple val(meta), path("*.gbf")    , emit: gbf_file
    tuple val(meta), path("19_table2asn.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_out = "${meta.mt_assembly_prefix}.sqn"
    def topology_mod = (circular?.toString()?.trim() == 'true') ? '[topology=circular] [completeness=complete]' : '[topology=linear]'
    def mgcode = meta.genetic_code ?: 2
    def effective_args = "-indir . -euk -J -t ${sample_sbt} -i ${sample_fa} -f ${sample_tbl} -w ${sample_cmt} -src-file ${sample_src} -o ${sample_out} -M n -j '[mgcode=${mgcode}] [location=mitochondrion] ${topology_mod}' -V vb -Z -W"

    """
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
        -W

    # Normalise optional files to stable required names for the Python parser.
    val_file=""
    for f in *.val; do [ -e "\$f" ] && { val_file="\$f"; break; }; done
    dr_file=""
    for f in *.dr; do [ -e "\$f" ] && { dr_file="\$f"; break; }; done
    [ -n "\${val_file}" ] && [ "\${val_file}" != "${meta.mt_assembly_prefix}.val" ] && cp "\${val_file}" "${meta.mt_assembly_prefix}.val"
    [ -e "${meta.mt_assembly_prefix}.val" ] || touch "${meta.mt_assembly_prefix}.val"
    [ -n "\${dr_file}" ] && [ "\${dr_file}" != "${meta.mt_assembly_prefix}.dr" ] && cp "\${dr_file}" "${meta.mt_assembly_prefix}.dr"
    [ -e "${meta.mt_assembly_prefix}.dr" ] || touch "${meta.mt_assembly_prefix}.dr"

    printf '%s\n' '<tr><td>Table2ASN</td><td><samp>${effective_args}</samp></td><td>Generates GenBank files and non-fatal per-sample validation reports for ${meta.id}.</td></tr>' > 19_table2asn.tool_params_mqcrow.html
    table2asn_version=\$(table2asn -version 2>&1 | head -1 | sed 's/.* //')
    printf '"%s":\n    table2asn: "%s"\n' "${task.process}" "\${table2asn_version}" > versions.yml
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
    : > \${prefix}.stats
    : > \${prefix}.dr
    : > \${prefix}.gbf
    printf '%s\n' '<tr><td>Table2ASN</td><td><samp>${effective_args}</samp></td><td>Generates GenBank files and non-fatal per-sample validation reports for ${meta.id}.</td></tr>' > 19_table2asn.tool_params_mqcrow.html
    printf '"%s":\n    table2asn: "stub"\n' "${task.process}" > versions.yml
    """
}
