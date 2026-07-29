// Fetch SEVERAL candidate mitogenome references for a sample instead of the first
// one findMitoReference happens to hit.
//
// MITOHIFI_FINDMITOREFERENCE walks the sample's NCBI lineage (species -> genus ->
// family -> ...) and stops at the first complete mitogenome it finds. For a taxon
// with no congeneric record that first hit is an arbitrary member of whatever rank
// the walk reached, and MitoHiFi then recruits reads against it -- a reference the
// sample's most divergent gene blocks map to poorly, which is how a clean-looking
// but gene-incomplete collapse happens (OG2102, Rouleina attrita, 7/13 PCGs).
//
// findMitoReference.py already supports -n (report N genomes, in order of
// identification), so this module is the same script asked for more of them; the
// choice between them is made by REFERENCE_RANK, on read evidence rather than on
// the order NCBI happened to return. Emits every candidate; always exits 0 with an
// empty candidate set on failure, so re-selection degrades to the single-reference
// behaviour instead of breaking the run.
process REFERENCE_CANDIDATES {
    tag "${meta.id} ${meta.reference_species_id ?: meta.nominal_species_id}"
    label 'process_single'
    label 'error_retry'
    secret secrets.NCBI_API_KEY ? "NCBI_API_KEY" : ""

    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    val(meta)

    output:
    tuple val(meta), path("candidates/*.fasta"), path("candidates/*.gb"), emit: candidates, optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix ?: meta.id}.reference_candidates_status.tsv"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }
    def args = task.ext.args ?: ''
    def ncbi_api_key = secrets.NCBI_API_KEY ? "--ncbi-api-key \$NCBI_API_KEY" : ""
    def n_candidates = (params.n_reference_candidates ?: 5) as int
    // Same species-name resolution as MITOHIFI_FINDMITOREFERENCE: prefer the
    // NCBI-valid name resolved at samplesheet creation, fall back to the raw
    // nominal name, and strip a trailing punctuation-only token that NCBI rejects.
    def ref_species = (meta.reference_species_id instanceof List)
        ? meta.reference_species_id.join('').trim()
        : (meta.reference_species_id ?: '').toString().trim()
    def query_species = (ref_species ?: meta.nominal_species_id).toString().replaceAll(/\s+/, ' ').replaceAll(/\s+[^A-Za-z]+$/, '').trim()
    def prefix = meta.mt_assembly_prefix ?: meta.id
    """
    export MPLBACKEND=Agg
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    mkdir -p .mplconfig candidates
    export NCBI_DELAY="\${NCBI_DELAY:-0.35}"

    set +e
    findMitoReference.py \\
        ${ncbi_api_key} \\
        --species "${query_species}" \\
        --outfolder candidates \\
        -n ${n_candidates} \\
        $args
    candidates_exit=\$?
    set -e

    n_found=\$(find candidates -maxdepth 1 -name '*.gb' -type f -size +0c | wc -l)
    if [ "\$candidates_exit" -ne 0 ] && [ '${task.attempt}' -le 3 ]; then
        # Preserve technical failures for the first three attempts so the error
        # strategy can retry them, matching MITOHIFI_FINDMITOREFERENCE.
        exit "\$candidates_exit"
    fi

    if [ "\$n_found" -eq 0 ]; then
        # No candidates is not fatal: the caller keeps the reference findMitoReference
        # already resolved. Drop any half-written pair so the optional output stays empty.
        rm -f candidates/*.fasta candidates/*.gb
        status=no_candidates
    else
        status=found
    fi
    printf 'sample\\tstatus\\tn_candidates\\texit_code\\tquery_species\\n%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        '${meta.id}' "\$status" "\$n_found" "\$candidates_exit" '${query_species}' \\
        > ${prefix}.reference_candidates_status.tsv

    printf '"%s":\\n    mitohifi: 3.2.3\\n' '${task.process}' > versions.yml
    """

    stub:
    def prefix = meta.mt_assembly_prefix ?: meta.id
    """
    mkdir -p candidates
    touch candidates/STUB001.fasta candidates/STUB001.gb
    printf 'sample\\tstatus\\tn_candidates\\texit_code\\tquery_species\\n${meta.id}\\tfound\\t1\\t0\\tstub\\n' > ${prefix}.reference_candidates_status.tsv
    printf '"%s":\\n    mitohifi: 3.2.3\\n' '${task.process}' > versions.yml
    """
}
