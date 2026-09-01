// Recover EMMA-dropped tRNAs in place from a tRNAscan-SE scan.
//
// EMMA's covariance model sometimes misses a tRNA that is physically present in
// the assembly. The gate (TRNA_RESCUE_GATE) only routes an assembly here when
// every missing REF gene is a tRNA and the whole 13-PCG + 2-rRNA core is present
// and ordered, so each target's insertion gap is well defined. TRNA_SCAN has
// already run tRNAscan-SE; rescue_trna.py parses that output and splices back the
// single best hit per target that passes an isotype/anticodon, score, length,
// gap-placement and zero-overlap guard, writing gene + tRNA lines into the .gff
// and .tbl exactly as EMMA writes its own tRNAs.
//
// Emits the SAME `results` bundle shape as EMMA so the annotation subworkflow
// mixes the rescued (FIX) and untouched (PASS) assemblies without rewiring.
//
// Fail-safe: rescue_trna.py guards every edit and always exits 0. A target that
// fails a guard (or an empty scan) is left untouched and the original bundle is
// re-emitted as-is, so the assembly still flows to QC and is held there exactly
// as before (unless the residual shortfall is within annotation_trna_tolerance).
process TRNA_RESCUE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    // bundle  = annotation bundle for a FIX assembly (EMMA, ND4L/ATP8-rescued)
    // scan    = TRNA_SCAN.out.tsv (tRNAscan-SE tabular output; may be empty)
    // targets = comma-list of REF tRNA names from the gate ('TP', 'TW,TA,TN', ...)
    tuple val(meta), path(bundle, stageAs: 'emma_in/*'), path(scan), val(targets)

    output:
    tuple val(meta), path("annotation/*"), emit: results
    tuple val(meta), path("annotation/mtdna_rescue/${meta.mt_assembly_prefix}.trna_rescue.status.txt"), emit: status
    tuple val(meta), path("09_trna_rescue.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
    def effective_args = "rescue_trna.py --annotation-dir annotation --targets ${targets} --scan-out ${scan} ${args}".replaceAll(/ +/, ' ').trim()
    """
    # Rebuild the annotation/ dir from the staged bundle (staged flat under
    # emma_in/ to keep it clear of this task's outputs).
    mkdir -p annotation
    cp -rL emma_in/* annotation/

    mkdir -p annotation/mtdna_rescue
    status_file=annotation/mtdna_rescue/${prefix}.trna_rescue.status.txt

    # rescue_trna.py guards every edit and exits 0 by design. The `|| echo` is a
    # second layer: an unexpected crash must still leave the bundle (already
    # copied above) intact and emitted, so the assembly keeps flowing to QC
    # exactly as if the rescue had never run.
    rescue_trna.py \\
        --annotation-dir annotation \\
        --targets ${targets} \\
        --scan-out ${scan} \\
        --status "\$status_file" \\
        ${args} || printf 'SKIP\\t-\\trescue_trna.py crashed\\n' > "\$status_file"

    [ -s "\$status_file" ] || printf 'SKIP\\t-\\tno status written\\n' > "\$status_file"

    status=\$(cut -f1 "\$status_file" | paste -sd, -)
    cat <<-END_TOOL_PARAMS > 09_trna_rescue.tool_params_mqcrow.html
    <tr><td>tRNA rescue</td><td><samp>${effective_args}</samp></td><td>Splices EMMA-dropped tRNAs (${targets}) for ${meta.id} from the tRNAscan-SE scan (outcome: \${status}); guarded by isotype/anticodon match, Infernal score, gap placement and a zero-overlap check.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
    """
    mkdir -p annotation/mtdna_rescue
    cp -rL emma_in/* annotation/ 2>/dev/null || true
    printf 'SKIP\\t-\\tstub\\n' > annotation/mtdna_rescue/${prefix}.trna_rescue.status.txt
    : > 09_trna_rescue.tool_params_mqcrow.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
