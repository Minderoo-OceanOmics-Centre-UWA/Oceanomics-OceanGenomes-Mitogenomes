// Recover an EMMA-dropped protein-coding gene (ND4L and/or ATP8) in place.
//
// EMMA finds these short genes but its rationalise_matches! overlap filter
// discards them against their longer neighbour (ND4, ATP6). The gene is present
// in the assembly, so this rebuilds the feature from the flanking-gene
// coordinates EMMA already produced: define the intergenic window, tblastn a
// small reference-protein set to fix the frame and guard identity/coverage,
// refine to a clean ORF, and splice matching gene/mRNA/CDS lines into the .gff
// and .tbl plus the per-gene cds/ and proteins/ FASTAs.
//
// Emits the SAME `results` bundle shape as EMMA so the annotation subworkflow
// mixes the rescued (FIX) and untouched (PASS) assemblies without rewiring.
//
// Fail-safe: rescue_emma_pcg.py guards every edit and always exits 0. A target
// that fails a guard is left untouched and the original EMMA bundle is re-emitted
// as-is, so the assembly still flows to LCA / species validation and is held at
// the QC gate exactly as before. Runs in the MITOS2 BioContainer (biopython +
// BLAST+, already used by the coral fixer).
process EMMA_GENE_RESCUE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // bundle  = EMMA.out.results for a FIX assembly
    // targets = comma-list from the gate ('ND4L', 'ATP8', 'ND4L,ATP8')
    // ref_faa = assets/rescue_pcg_refs.faa
    tuple val(meta), path(bundle, stageAs: 'emma_in/*'), val(targets), path(ref_faa)

    output:
    tuple val(meta), path("annotation/*"), emit: results
    tuple val(meta), path("annotation/mtdna_rescue/${meta.mt_assembly_prefix}.rescue.status.txt"), emit: status
    tuple val(meta), path("08_emma_rescue.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def code   = task.ext.code ?: meta.genetic_code
    def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
    def effective_args = "rescue_emma_pcg.py --annotation-dir annotation --targets ${targets} --ref-faa ${ref_faa} --code ${code} ${args}".replaceAll(/ +/, ' ').trim()
    """
    # Rebuild the annotation/ dir from the staged EMMA bundle (staged flat under
    # emma_in/ to keep it clear of this task's outputs).
    mkdir -p annotation
    cp -rL emma_in/* annotation/

    mkdir -p annotation/mtdna_rescue
    status_file=annotation/mtdna_rescue/${prefix}.rescue.status.txt

    # rescue_emma_pcg.py guards every edit and exits 0 by design. The `|| true`
    # is a second layer: an unexpected crash must still leave the EMMA bundle
    # (already copied above) intact and emitted, so the assembly keeps flowing to
    # LCA / QC exactly as if the rescue had never run.
    rescue_emma_pcg.py \\
        --annotation-dir annotation \\
        --targets ${targets} \\
        --ref-faa ${ref_faa} \\
        --code ${code} \\
        --status "\$status_file" \\
        ${args} || printf 'SKIP\\t-\\trescue_emma_pcg.py crashed\\n' > "\$status_file"

    [ -s "\$status_file" ] || printf 'SKIP\\t-\\tno status written\\n' > "\$status_file"

    status=\$(cut -f1 "\$status_file" | paste -sd, -)
    cat <<-END_TOOL_PARAMS > 08_emma_rescue.tool_params_mqcrow.html
    <tr><td>EMMA gene rescue</td><td><samp>${effective_args}</samp></td><td>Recovers EMMA-dropped ${targets} for ${meta.id} from the assembly (outcome: \${status}); guarded by reference BLAST identity/coverage and an ORF check.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/^Python //')
        tblastn: \$(tblastn -version 2>/dev/null | sed -n 's/^tblastn: //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
    """
    mkdir -p annotation/cds annotation/proteins annotation/mtdna_rescue
    cp -rL emma_in/* annotation/ 2>/dev/null || true
    printf 'SKIP\\t-\\tstub\\n' > annotation/mtdna_rescue/${prefix}.rescue.status.txt
    : > 08_emma_rescue.tool_params_mqcrow.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
        tblastn: "2.17.0"
    END_VERSIONS
    """
}
