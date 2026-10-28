// Extend EMMA protein-coding genes whose stop codon EMMA cut one codon short.
//
// Where EMMA places a PCG stop depends on where the assembler linearised the
// circle: the same genome at two rotations can come back with CO1 ending on the
// encoded TAA or one codon earlier, on a sense codon EMMA labels a "putative stop
// completed by 3' A residues". NOVA_260925_LA held 11 Glaucosoma hebraicum on
// exactly this (9 CO1, 2 ND1); four were byte-identical genomes to passing
// batch-mates, only rotated. The vertebrate path had no ORF repair at all.
//
// emma_stop_repair.py moves a CDS's 3' end downstream by at most 3 nt, onto a
// TAA/TAG the assembly itself encodes, and only when the added bases sit inside
// no other annotated feature (which is what leaves EMMA's genuine poly-A stops
// alone). It rewrites the .gff, .tbl and that gene's cds/ + proteins/ FASTAs.
//
// Runs on EVERY EMMA bundle (it is cheap and decides per gene); a bundle with
// nothing to repair is re-emitted byte for byte. Emits the SAME `results` shape
// as EMMA so the subworkflow swaps it in without rewiring, plus CO1 so BLAST/LCA
// sees the repaired sequence.
//
// Fail-safe: the script guards every edit and always exits 0. A crash still
// leaves the bundle (copied first) intact and emitted.
process EMMA_STOP_REPAIR {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    // bundle = annotation bundle after the ND4L/ATP8 and tRNA rescues
    tuple val(meta), path(bundle, stageAs: 'emma_in/*')

    output:
    tuple val(meta), path("annotation/*"), emit: results
    tuple val(meta), path("annotation/cds/*CO1*.fa"), emit: co1_sequences, optional: true
    tuple val(meta), path("annotation/mtdna_rescue/${meta.mt_assembly_prefix}.stop_repair.status.txt"), emit: status
    tuple val(meta), path("09b_emma_stop_repair.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def code   = task.ext.code ?: meta.genetic_code
    def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
    def effective_args = "emma_stop_repair.py --annotation-dir annotation --code ${code} ${args}".replaceAll(/ +/, ' ').trim()
    """
    # emma_stop_repair.py imports orf_utils.py -- named here on purpose so that
    # Nextflow hashes it (it only hashes bin/ scripts that appear in the command).
    mkdir -p annotation
    cp -rL emma_in/* annotation/

    mkdir -p annotation/mtdna_rescue
    status_file=annotation/mtdna_rescue/${prefix}.stop_repair.status.txt

    emma_stop_repair.py \\
        --annotation-dir annotation \\
        --code ${code} \\
        --status "\$status_file" \\
        ${args} || printf 'SKIP\\t-\\temma_stop_repair.py crashed\\n' > "\$status_file"

    [ -s "\$status_file" ] || printf 'SKIP\\t-\\tno status written\\n' > "\$status_file"

    outcome=\$(awk -F'\\t' '{print \$1 (\$2 != "-" ? " " \$2 : "")}' "\$status_file" | paste -sd, -)
    cat <<-END_TOOL_PARAMS > 09b_emma_stop_repair.tool_params_mqcrow.html
    <tr><td>EMMA stop repair</td><td><samp>${effective_args}</samp></td><td>Extends EMMA PCGs cut short of an encoded TAA/TAG for ${meta.id} (outcome: \${outcome}); at most 3 nt, never into another feature.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
    """
    mkdir -p annotation/cds annotation/proteins annotation/mtdna_rescue
    cp -rL emma_in/* annotation/ 2>/dev/null || true
    printf 'NONE\\t-\\tstub\\n' > annotation/mtdna_rescue/${prefix}.stop_repair.status.txt
    : > 09b_emma_stop_repair.tool_params_mqcrow.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
