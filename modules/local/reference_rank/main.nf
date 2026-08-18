// Choose the best reference for a sample from the REFERENCE_CANDIDATES set, by
// mapping a subsample of the sample's own reads against each candidate.
//
// MitoHiFi builds its assembly by recruiting reads that map to the reference, so
// "how many of this sample's reads does this candidate recruit" is exactly the
// property that decides whether the assembly succeeds -- and unlike the finished
// assembly, it is measurable *before* assembly, from reads alone. That is what
// makes re-selection possible at the point where it can still change the outcome.
//
// Emits the winning fasta/gb plus a ranking TSV showing every candidate and its
// score, so a curator can see why one was chosen. Runs in the MitoHiFi container,
// which already ships minimap2. Always exits 0.
//
// The chosen_reference outputs are EMPTY when the script declined to substitute --
// no candidate recruited any reads, no reads were subsampled, or there was a single
// unmapped candidate. Callers must test them for size and fall back to the reference
// findMitoReference resolved; an empty file here is a verdict, not a failure. See the
// header of bin/rank_reference_candidates.py for why a zero-scoring candidate set is
// not a tie to be broken.
process REFERENCE_RANK {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    // candidate_fastas / candidate_gbs are the full candidate set from
    // REFERENCE_CANDIDATES; reads are the sample's (already concatenated) reads.
    tuple val(meta), path(reads), path(candidate_fastas), path(candidate_gbs)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.chosen_reference.fasta"),
                     path("${meta.mt_assembly_prefix}.chosen_reference.gb"), emit: reference
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference_ranking.tsv"), emit: ranking
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // HiFi reads need map-hifi; Illumina/HiC are short reads and need sr. Getting
    // this wrong would score every candidate near zero and make the ranking noise.
    def preset = meta.sequencing_type == 'hifi' ? 'map-hifi' : 'sr'
    def subsample = (params.reference_rank_read_subsample ?: 50000) as int
    """
    mkdir -p candidates
    # Candidates are staged flat; the script pairs <acc>.fasta with <acc>.gb by stem.
    for f in ${candidate_fastas} ${candidate_gbs}; do
        cp -L "\$f" candidates/ 2>/dev/null || true
    done

    rank_reference_candidates.py \\
        --candidate-dir candidates \\
        --reads ${reads} \\
        --prefix ${meta.mt_assembly_prefix} \\
        --preset ${preset} \\
        --subsample-reads ${subsample} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>/dev/null)
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    // Non-empty on purpose: empty chosen_reference files now MEAN "declined to
    // substitute", and this stub's ranking says chosen=yes. touch-ing them would
    // make every stub run silently take the fallback path it claims not to.
    """
    printf '>STUB001\\nACGT\\n' > ${meta.mt_assembly_prefix}.chosen_reference.fasta
    printf 'LOCUS       STUB001                 4 bp    DNA     circular UNK\\n//\\n' > ${meta.mt_assembly_prefix}.chosen_reference.gb
    printf 'accession\\tscore\\tmatched_bases\\tmapped_reads\\tchosen\\tstatus\\nSTUB001\\t1.0\\t100\\t10\\tyes\\tscored\\n' > ${meta.mt_assembly_prefix}.reference_ranking.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "2.24-r1122"
        python: "3.10.0"
    END_VERSIONS
    """
}
