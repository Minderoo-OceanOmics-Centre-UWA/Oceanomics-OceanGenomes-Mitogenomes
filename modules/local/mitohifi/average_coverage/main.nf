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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    bam="${coverage_mapping}/HiFi-vs-final_mitogenome.sorted.bam"

    if [ -s "\$bam" ]; then
        samtools depth -aa "\$bam" | awk -v sample="${prefix}" 'BEGIN {
            OFS="\\t"
            print "sample", "contig_id", "length", "covered_bases", "avg_coverage"
        }
        {
            contig=\$1
            sum += \$3
            n += 1
            if (\$3 > 0) covered += 1
        }
        END {
            if (n > 0) {
                printf "%s\\t%s\\t%d\\t%d\\t%.6f\\n", sample, contig, n, covered, sum / n
            } else {
                printf "%s\\tNA\\t0\\t0\\tNA\\n", sample
            }
        }' > ${prefix}.coverage.tsv
    else
        printf "sample\\tcontig_id\\tlength\\tcovered_bases\\tavg_coverage\\n%s\\tNA\\t0\\t0\\tNA\\n" "${prefix}" > ${prefix}.coverage.tsv
    fi

    awk 'BEGIN { FS=OFS="\\t" }
        FNR == NR {
            if (FNR > 1 && \$2 != "NA") {
                cov[\$2] = \$5
                contig = \$2
                suffix = "_rotated"
                if (length(contig) > length(suffix) && substr(contig, length(contig) - length(suffix) + 1) == suffix) {
                    contig = substr(contig, 1, length(contig) - length(suffix))
                }
                cov[contig] = \$5
                if (default_cov == "") default_cov = \$5
            }
            next
        }
        /^#/ {
            print
            next
        }
        !header_seen {
            print \$0, "avg_coverage"
            header_seen = 1
            next
        }
        {
            c = (\$1 in cov) ? cov[\$1] : (\$1 == "final_mitogenome" ? default_cov : "")
            print \$0, c
        }' ${prefix}.coverage.tsv ${contigs_stats} > ${prefix}.contigs_stats.with_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    cp ${contigs_stats} ${prefix}.contigs_stats.with_coverage.tsv
    printf "sample\\tcontig_id\\tlength\\tcovered_bases\\tavg_coverage\\n%s\\tNA\\t0\\t0\\tNA\\n" "${prefix}" > ${prefix}.coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "1.20"
    END_VERSIONS
    """
}
