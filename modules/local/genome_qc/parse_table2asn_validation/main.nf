process PARSE_TABLE2ASN_VALIDATION {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0'

    input:
    tuple val(meta), path(val_file), path(discrepancy_file), val(circular)

    output:
    tuple val(meta), path("${prefix}.qc_flags.tsv"), emit: qc_flags
    tuple val(meta), path("${prefix}.table2asn_findings.tsv"), emit: findings
    tuple val(meta), path("${prefix}.table2asn_status.tsv"), emit: status
    path "versions.yml", emit: versions

    script:
    prefix = meta.mt_assembly_prefix ?: meta.id
    """
    parse_table2asn_validation.py \\
        --sample '${prefix}' \\
        --circular '${circular}' \\
        --val ${val_file} \\
        --dr ${discrepancy_file} \\
        --findings ${prefix}.table2asn_findings.tsv \\
        --status ${prefix}.table2asn_status.tsv \\
        --qc-flags ${prefix}.qc_flags.tsv
    printf '"%s":\n    python: "%s"\n' '${task.process}' "\$(python --version | sed 's/Python //')" > versions.yml
    """

    stub:
    prefix = meta.mt_assembly_prefix ?: meta.id
    """
    printf 'sample\tsource\tseverity\tcode\tmessage\n' > ${prefix}.table2asn_findings.tsv
    printf 'sample\tstatus\treject_count\terror_count\twarning_count\tinfo_count\tfatal_discrepancy_count\tnostop_count\tblocking_codes\twarning_codes\n${prefix}\tPASS\t0\t0\t0\t0\t0\t0\t\t\n' > ${prefix}.table2asn_status.tsv
    printf 'sample\tcircular\tstatus\treject_count\terror_count\twarning_count\tfatal_discrepancy_count\tnostop_count\tnostop_features\n${prefix}\t${circular}\tPASS\t0\t0\t0\t0\t0\t\n' > ${prefix}.qc_flags.tsv
    printf '"%s":\n    python: "stub"\n' '${task.process}' > versions.yml
    """
}
