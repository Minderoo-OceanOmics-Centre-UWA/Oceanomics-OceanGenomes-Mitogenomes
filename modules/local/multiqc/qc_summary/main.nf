process QC_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(species_file), path(proceed_file), path(annotation_stats_csv)

    output:
    tuple val(meta), path("${meta.id}.qc_summary.tsv"), emit: table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    species=\$(tr -d '\r' < ${species_file} | tr -d '\n')
    proceed=\$(tr -d '\r' < ${proceed_file} | tr -d '\n')
    # Extract columns from annotation CSV if present
    ann_passed=NA
    missing=NA
    if [ -s "${annotation_stats_csv}" ]; then
      # Find column indices
      pass_col=\$(awk -F"," 'NR==1{for(i=1;i<=NF;i++){if(\$i=="passed"){print i; exit}}}' ${annotation_stats_csv})
      miss_col=\$(awk -F"," 'NR==1{for(i=1;i<=NF;i++){if(\$i=="missing_genes"||\$i=="num_missing"){print i; exit}}}' ${annotation_stats_csv})
      if [ -n "\${pass_col}" ]; then
        ann_passed=\$(awk -F"," -v c=\${pass_col} 'NR==2{print \$c}' ${annotation_stats_csv})
      fi
      if [ -n "\${miss_col}" ]; then
        missing=\$(awk -F"," -v c=\${miss_col} 'NR==2{print \$c}' ${annotation_stats_csv})
      fi
    fi
    out=${meta.id}.qc_summary.tsv
    {
      echo -e "sample\\tspecies\\tproceed_qc\\tannotation_passed\\tmissing_genes"
      echo -e "${meta.id}\\t\${species}\\t\${proceed}\\t\${ann_passed}\\t\${missing}"
    } > \${out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qc_summary: "1.1.0"
    END_VERSIONS
    """

    stub:
    """
    out=${meta.id}.qc_summary.tsv
    {
      echo -e "sample\\tspecies\\tproceed_qc\\tannotation_passed\\tmissing_genes"
      echo -e "${meta.id}\\ttest_species\\ttrue\\tyes\\t0"
    } > \${out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qc_summary: "stub"
    END_VERSIONS
    """
}
