process OATK_CHECK {
    tag "$meta.id"
    label 'process_low'

    // Reuse the pinned MitoHiFi image (blastn + python): the same FASTA-vs-reference
    // circularity/anomaly engine the GetOrganelle check uses. Oatk produces no
    // read-mapping BAM, so the BAM-based MitoHiFi check does not apply here; this
    // reference-coverage + self-align check does, and needs no extra container.
    container "${params.mitohifi_container}"

    input:
    // fasta  : the assembled oatk contig (<prefix>.fasta).
    // gfa    : oatk's assembly graph. A link line joining a segment to itself
    //          (L <seg> + <seg> +) is a closed loop -> circular; this is oatk's own
    //          circularity signal, fed to the check as --getorg-circular.
    // ref_gb : related-species reference GenBank, or assets/placeholders/NO_REFERENCE.gb when the
    //          sample reached oatk via the no-reference path (check treats a zero-length
    //          reference as absent, so the length/repeat/concatemer screen still runs).
    tuple val(meta), path(fasta), path(gfa), path(reference_gb)

    output:
    // Emit the GetOrganelle-check schema AND name (<prefix>.getorg_check.tsv) so the
    // assembly summary applies its getorg_circular_override (final_verdict_circular ->
    // circularised) and reads anomaly_type/length_anomaly by column name, giving the
    // oatk contig the SAME summary treatment as a GetOrganelle assembly. The parent
    // workflow reads final_verdict_circular / anomaly_type by name, so the filename is
    // immaterial there.
    tuple val(meta), path("${prefix}.getorg_check.tsv")            , emit: evidence
    tuple val(meta), path("09_oatk_check.tool_params_mqcrow.html") , emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    # Oatk circularity: a GFA link line joining a segment to itself (L <seg> + <seg> +)
    # closes the contig into a circle. No such self-link -> treat as non-circular; the
    # reference-coverage test below can still upgrade it, mirroring the GetOrganelle path.
    circ=false
    if [ -s ${gfa} ] && awk -F'\\t' '\$1=="L" && \$2==\$4 {found=1} END{exit !found}' ${gfa}; then
        circ=true
    fi

    check_getorganelle.py \\
        --fasta ${fasta} \\
        --reference-gb ${reference_gb} \\
        --sample ${prefix} \\
        --getorg-circular \$circ \\
        --threads ${task.cpus} \\
        --out-evidence ${prefix}.getorg_check.tsv \\
        ${args}

    cat <<-END_TOOL_PARAMS > 09_oatk_check.tool_params_mqcrow.html
    <tr><td>Oatk circularity check</td><td><samp>check_getorganelle.py --reference-gb (contig-vs-reference coverage) + self-align length/repeat screen; circularity seeded from the oatk GFA self-link</samp></td><td>Runs the same FASTA-based circularity/anomaly check on the reference-free oatk fallback contig for ${meta.id}, so it flows into annotation with a real topology verdict like the MitoHiFi / GetOrganelle assemblies.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version 2>/dev/null | head -n1 | sed 's/^blastn: //')
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    printf "sample\\tgetorg_circular\\tnum_records\\treference_coverage\\treference_max_gap\\tcircular_by_reference\\tfinal_verdict_circular\\tcircular_corrected\\tassembly_length\\treference_length\\tlength_ratio\\texcess_bp\\tlength_anomaly\\ttandem_repeat\\trepeat_region\\tanomaly_type\\tsuggested_trim_region\\tcuration_suggestion\\tnote\\n%s\\ttrue\\t1\\tNA\\tNA\\tNA\\ttrue\\tno\\tNA\\tNA\\tNA\\tNA\\tNA\\tno\\tNA\\tnone\\tNA\\tnone\\tstub\\n" "${prefix}" > ${prefix}.getorg_check.tsv

    cat <<-END_TOOL_PARAMS > 09_oatk_check.tool_params_mqcrow.html
    <tr><td>Oatk circularity check</td><td><samp>check_getorganelle.py (stub)</samp></td><td>Oatk circularity check for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: "stub"
        python: "stub"
    END_VERSIONS
    """
}
