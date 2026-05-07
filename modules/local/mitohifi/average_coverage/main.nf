process MITOHIFI_AVERAGE_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
        'quay.io/biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(contigs_stats), path(coverage_mapping)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.contigs_stats.with_coverage.tsv"), emit: stats
    tuple val(meta), path("${meta.mt_assembly_prefix}.coverage.tsv"), emit: coverage
    tuple val(meta), path("06_mitohifi_average_coverage.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    def effective_args = "samtools depth -aa ${coverage_mapping}/HiFi-vs-final_mitogenome.sorted.bam"
    """
    bam="${coverage_mapping}/HiFi-vs-final_mitogenome.sorted.bam"

    if [ -s "\$bam" ]; then
        samtools depth -aa "\$bam" | awk -v sample="${prefix}" 'BEGIN {
            OFS="\\t"
            print "sample", "contig_id", "length", "covered_bases", "avg_coverage", "coverage_cv"
        }
        {
            contig=\$1
            depth=\$3
            sum += depth
            sumsq += depth * depth
            n += 1
            if (depth > 0) covered += 1
        }
        END {
            if (n > 0) {
                mean = sum / n
                variance = (sumsq / n) - (mean * mean)
                if (variance < 0) variance = 0
                cv = mean > 0 ? sqrt(variance) / mean : "NA"
                printf "%s\\t%s\\t%d\\t%d\\t%.6f\\t%s\\n", sample, contig, n, covered, mean, cv
            } else {
                printf "%s\\tNA\\t0\\t0\\tNA\\tNA\\n", sample
            }
        }' > ${prefix}.coverage.tsv
    else
        printf "sample\\tcontig_id\\tlength\\tcovered_bases\\tavg_coverage\\tcoverage_cv\\n%s\\tNA\\t0\\t0\\tNA\\tNA\\n" "${prefix}" > ${prefix}.coverage.tsv
    fi

    awk 'BEGIN { FS=OFS="\\t" }
        FNR == NR {
            if (FNR > 1 && \$2 != "NA") {
                cov[\$2] = \$5
                cv[\$2] = \$6
                contig = \$2
                suffix = "_rotated"
                if (length(contig) > length(suffix) && substr(contig, length(contig) - length(suffix) + 1) == suffix) {
                    contig = substr(contig, 1, length(contig) - length(suffix))
                }
                cov[contig] = \$5
                cv[contig] = \$6
                if (default_cov == "") default_cov = \$5
                if (default_cv == "") default_cv = \$6
            }
            next
        }
        /^#/ {
            print
            next
        }
        !header_seen {
            print \$0, "avg_coverage", "coverage_cv"
            header_seen = 1
            next
        }
        {
            c = (\$1 in cov) ? cov[\$1] : (\$1 == "final_mitogenome" ? default_cov : "")
            v = (\$1 in cv) ? cv[\$1] : (\$1 == "final_mitogenome" ? default_cv : "")
            print \$0, c, v
        }' ${prefix}.coverage.tsv ${contigs_stats} > ${prefix}.contigs_stats.with_coverage.tsv

    cat <<-END_TOOL_PARAMS > 06_mitohifi_average_coverage.tool_params_mqcrow.html
    <tr><td>MitoHiFi Average Coverage</td><td><samp>${effective_args}</samp></td><td>Adds average coverage from the MitoHiFi coverage BAM to contig statistics for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | awk 'NR==1 {print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    def effective_args = "samtools depth -aa ${coverage_mapping}/HiFi-vs-final_mitogenome.sorted.bam"
    """
    cp ${contigs_stats} ${prefix}.contigs_stats.with_coverage.tsv
    printf "sample\\tcontig_id\\tlength\\tcovered_bases\\tavg_coverage\\tcoverage_cv\\n%s\\tNA\\t0\\t0\\tNA\\tNA\\n" "${prefix}" > ${prefix}.coverage.tsv

    cat <<-END_TOOL_PARAMS > 06_mitohifi_average_coverage.tool_params_mqcrow.html
    <tr><td>MitoHiFi Average Coverage</td><td><samp>${effective_args}</samp></td><td>Adds average coverage from the MitoHiFi coverage BAM to contig statistics for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "1.20"
    END_VERSIONS
    """
}
