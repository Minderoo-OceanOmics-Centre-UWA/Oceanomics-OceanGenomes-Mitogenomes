process MITOHIFI_CHECK_CIRCULARITY {
    tag "$meta.id"
    label 'process_low'

    // Reuse the pinned MitoHiFi image: it ships the exact minimap2 + samtools +
    // python the read-spanning junction test needs, so the check runs the same
    // aligner MitoHiFi used and needs no extra container.
    container "${params.mitohifi_container}"

    input:
    tuple val(meta), path(contigs_stats), path(fasta), path(coverage_mapping), path(gb)

    output:
    // Corrected contig-stats: drop-in replacement for the average-coverage table
    // (same basename) with was_circular forced True when the re-test proves the
    // assembly circular. Written under checked/ to avoid clashing with the staged
    // input of the same name.
    tuple val(meta), path("checked/${meta.mt_assembly_prefix}.contigs_stats.with_coverage.tsv"), emit: stats
    tuple val(meta), path("${meta.mt_assembly_prefix}.circularity_check.tsv")                   , emit: evidence
    tuple val(meta), path("07_check_circularity.tool_params_mqcrow.html")                       , emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    def bam = "${coverage_mapping}/HiFi-vs-final_mitogenome.sorted.bam"
    """
    mkdir -p checked

    check_circularity.py \\
        --fasta ${fasta} \\
        --contigs-stats ${contigs_stats} \\
        --bam ${bam} \\
        --gb ${gb} \\
        --sample ${prefix} \\
        --threads ${task.cpus} \\
        --out-stats checked/${prefix}.contigs_stats.with_coverage.tsv \\
        --out-evidence ${prefix}.circularity_check.tsv \\
        ${args}

    cat <<-END_TOOL_PARAMS > 07_check_circularity.tool_params_mqcrow.html
    <tr><td>Circularity check</td><td><samp>check_circularity.py --bam HiFi-vs-final_mitogenome.sorted.bam (doubled-reference remap) + hifiasm c/l flag; self-align vs reference length for control-region repeats</samp></td><td>Re-tests MitoHiFi non-circular calls for ${meta.id} (junction-spanning reads + hifiasm flag) and corrects was_circular in place; also flags over-length assemblies whose excess is a control-region tandem repeat for manual curation.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>/dev/null)
        samtools: \$(samtools --version | awk 'NR==1 {print \$2}')
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    mkdir -p checked
    cp ${contigs_stats} checked/${prefix}.contigs_stats.with_coverage.tsv
    printf "sample\\tsource_contig\\tmitohifi_was_circular\\thifiasm_contig_circular\\tjunction_spanning_reads\\tmapped_reads\\tmin_spanning_reads\\tmin_overhang\\tread_spanning_circular\\tfinal_verdict_circular\\tflag_corrected\\tnote\\tassembly_length\\treference_length\\tlength_ratio\\texcess_bp\\tlength_anomaly\\ttandem_repeat\\trepeat_region\\trepeat_unit_bp\\trepeat_copies\\trepeat_identity\\trepeat_in_control_region\\tanomaly_type\\tsuggested_trim_region\\tcuration_suggestion\\n%s\\tNA\\tNA\\tNA\\tNA\\tNA\\t2\\t200\\tNA\\tNA\\tno\\tstub\\tNA\\tNA\\tNA\\tNA\\tNA\\tno\\tNA\\tNA\\tNA\\tNA\\tNA\\tnone\\tNA\\tnone\\n" "${prefix}" > ${prefix}.circularity_check.tsv

    cat <<-END_TOOL_PARAMS > 07_check_circularity.tool_params_mqcrow.html
    <tr><td>Circularity check</td><td><samp>check_circularity.py (stub)</samp></td><td>Re-tests MitoHiFi non-circular calls for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "stub"
        samtools: "stub"
        python: "stub"
    END_VERSIONS
    """
}
