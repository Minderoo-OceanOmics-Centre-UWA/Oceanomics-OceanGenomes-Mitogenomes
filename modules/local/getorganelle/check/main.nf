process GETORGANELLE_CHECK {
    tag "$meta.id"
    label 'process_low'

    // Reuse the pinned MitoHiFi image for blastn + python (same tooling as the
    // HiFi circularity check); no extra container needed.
    container "${params.mitohifi_container}"

    input:
    tuple val(meta), path(fasta), path(reference_gb)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.getorg_check.tsv")            , emit: evidence
    tuple val(meta), path("08_getorganelle_check.tool_params_mqcrow.html")          , emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    // GetOrganelle log verdict carried on meta.circular (true/false/null). The
    // script only relabels false/unknown scaffolds the reference confirms circular.
    def circ = (meta.circular == null) ? 'null' : meta.circular.toString()
    // An empty placeholder GenBank (assets/NO_REFERENCE.gb) is passed for samples
    // with no findMitoReference; the script treats a zero-length reference as absent.
    """
    check_getorganelle.py \\
        --fasta ${fasta} \\
        --reference-gb ${reference_gb} \\
        --sample ${prefix} \\
        --getorg-circular ${circ} \\
        --threads ${task.cpus} \\
        --out-evidence ${prefix}.getorg_check.tsv \\
        ${args}

    cat <<-END_TOOL_PARAMS > 08_getorganelle_check.tool_params_mqcrow.html
    <tr><td>GetOrganelle check</td><td><samp>check_getorganelle.py --reference-gb (scaffold-vs-reference coverage) + self-align length/repeat screen</samp></td><td>Re-tests non-circular GetOrganelle scaffolds for ${meta.id} against the related reference (full reference coverage = a complete circle linearised elsewhere) and screens length/tandem-repeat anomalies; corrects the circular verdict when confirmed.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version 2>/dev/null | head -n1 | sed 's/^blastn: //')
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    printf "sample\\tgetorg_circular\\tnum_records\\treference_coverage\\treference_max_gap\\tcircular_by_reference\\tfinal_verdict_circular\\tcircular_corrected\\tassembly_length\\treference_length\\tlength_ratio\\texcess_bp\\tlength_anomaly\\ttandem_repeat\\trepeat_region\\tanomaly_type\\tsuggested_trim_region\\tcuration_suggestion\\tnote\\n%s\\tNA\\t1\\tNA\\tNA\\tNA\\tNA\\tno\\tNA\\tNA\\tNA\\tNA\\tNA\\tno\\tNA\\tnone\\tNA\\tnone\\tstub\\n" "${prefix}" > ${prefix}.getorg_check.tsv

    cat <<-END_TOOL_PARAMS > 08_getorganelle_check.tool_params_mqcrow.html
    <tr><td>GetOrganelle check</td><td><samp>check_getorganelle.py (stub)</samp></td><td>Re-tests non-circular GetOrganelle scaffolds for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: "stub"
        python: "stub"
    END_VERSIONS
    """
}
