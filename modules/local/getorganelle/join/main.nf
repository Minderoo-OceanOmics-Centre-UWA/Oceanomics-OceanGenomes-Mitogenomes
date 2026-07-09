// Reference-guided join of a multi-scaffold GetOrganelle result.
//
// GetOrganelle sometimes returns >=2 scaffolds that together carry the full gene
// complement but that it could not disentangle into one circle (typically a control-
// region repeat the short reads cannot span). The reseed pass closes many of these;
// the survivors were, until now, only blindly concatenated by SANITISE_FASTA to
// scrape a species ID -- which also severs whichever gene straddles the seam and
// leaves minus-strand scaffolds reversed, costing MITOS2 annotation calls.
//
// This module runs bin/join_scaffolds_by_reference.py: it BLASTs each scaffold
// against the related reference (the same reference GETORGANELLE_CHECK consumes),
// orders the scaffolds by reference coordinate, orients them to the reference (+)
// strand, and joins them with N-spacers across genuine gaps. The single-record
// result then flows into GETORGANELLE_CHECK exactly as a one-scaffold result would.
//
// Only multi-scaffold samples WITH a real reference are routed here; the output is
// named "<assembly_prefix>_rgj" (reference_guided_join) so it is distinguishable
// from a first-pass/reseed assembly and from a SANITISE_FASTA "_concat" fallback.
process GETORGANELLE_JOIN {
    tag "$meta.id"
    label 'process_low'

    // Reuse the pinned MitoHiFi image for blastn + python (same tooling as
    // GETORGANELLE_CHECK / the HiFi circularity check); no extra container needed.
    container "${params.mitohifi_container}"

    input:
    tuple val(meta), path(fasta), path(reference_gb)

    output:
    // Named after the assembly evaluated (fasta basename) + "_rgj", not
    // meta.mt_assembly_prefix: the reseed pass carries its "reseed" suffix only on the
    // fasta filename, so keying off meta would drop it (getorg1770_rgj instead of
    // getorg1770reseed_rgj) and misfile the result under the first-pass prefix.
    tuple val(meta), path("*_rgj.{fa,fasta,fna}")        , emit: fasta
    // Pass the reference through so GETORGANELLE_CHECK can be re-paired with it
    // without re-consuming the (single-use, partial) reference channel.
    tuple val(meta), path(reference_gb)                  , emit: reference
    tuple val(meta), path("*.scaffold_join.tsv")         , emit: evidence
    tuple val(meta), path("09_getorganelle_join.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Strip only the FASTA extension (not simpleName, which strips every dot-field
    // and would collapse OG1946.ilmn.251215.getorg1770reseed down to OG1946).
    def base = fasta.name.replaceAll(/\.(fa|fasta|fna)$/, '')
    def prefix = "${base}_rgj"
    """
    join_scaffolds_by_reference.py \\
        --fasta ${fasta} \\
        --reference-gb ${reference_gb} \\
        --sample ${prefix} \\
        --out-fasta ${prefix}.fasta \\
        --out-evidence ${prefix}.scaffold_join.tsv \\
        ${args}

    cat <<-END_TOOL_PARAMS > 09_getorganelle_join.tool_params_mqcrow.html
    <tr><td>GetOrganelle scaffold join</td><td><samp>join_scaffolds_by_reference.py --reference-gb (dc-megablast scaffold ordering/orientation + N-gap join)</samp></td><td>Reference-guided order/orient/join of the multi-scaffold GetOrganelle result for ${meta.id} into a single best-effort molecule (flagged reference_guided_join).</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version 2>/dev/null | head -n1 | sed 's/^blastn: //')
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def base = fasta.name.replaceAll(/\.(fa|fasta|fna)$/, '')
    def prefix = "${base}_rgj"
    """
    printf ">%s\\nACGT\\n" "${prefix}" > ${prefix}.fasta
    printf "sample\\tmethod\\tn_input_scaffolds\\tn_placed\\tn_unplaced\\tscaffold_order\\tseams\\tjoined_length\\treference_length\\tlength_ratio\\treference_coverage\\torigin_wrap_gap_bp\\tlooks_circular\\tnote\\n%s\\treference_guided_join\\t2\\t2\\t0\\tNA\\tNA\\t4\\tNA\\tNA\\tNA\\tNA\\tNA\\tstub\\n" "${prefix}" > ${prefix}.scaffold_join.tsv

    cat <<-END_TOOL_PARAMS > 09_getorganelle_join.tool_params_mqcrow.html
    <tr><td>GetOrganelle scaffold join</td><td><samp>join_scaffolds_by_reference.py (stub)</samp></td><td>Reference-guided join for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: "stub"
        python: "stub"
    END_VERSIONS
    """
}
