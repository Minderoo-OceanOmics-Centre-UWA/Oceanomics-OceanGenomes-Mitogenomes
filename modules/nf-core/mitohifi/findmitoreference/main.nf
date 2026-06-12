process MITOHIFI_FINDMITOREFERENCE {
    tag "${meta.id} ${meta.reference_species_id ?: meta.nominal_species_id}"
    label 'process_single'
    label 'error_retry'
    secret secrets.NCBI_API_KEY ? "NCBI_API_KEY" : ""

    // NOTE: An optional NCBI API key raises the NCBI request-rate limit and
    // reduces transient lookup failures. Set it with Nextflow secrets:
    //   nextflow secrets set NCBI_API_KEY <key>
    // See https://www.nextflow.io/docs/latest/secrets.html for more information.

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    val(meta)

    output:
    tuple val(meta), path("*.fasta"), path("*.gb")  , emit: reference
    tuple val(meta), path("04_mitohifi_findmitoreference.tool_params_mqcrow.html"), emit: tool_params
    tuple val(meta), path("versions.yml")           , emit: versions_tuple
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }
    def args = task.ext.args ?: ''
    def ncbi_api_key = secrets.NCBI_API_KEY ? "--ncbi-api-key \$NCBI_API_KEY" : ""
    // Prefer the NCBI-valid name resolved at samplesheet creation; fall back to
    // the raw nominal name for older samplesheets that lack the column. The raw
    // nominal name is frequently rejected by NCBI ('No such species in NCBI!').
    // Older samplesheets lack the reference_species_id column; nf-schema fills the
    // missing meta value with an empty list (not null), so coerce to a clean string
    // before testing it and fall back to the raw nominal name when it's empty.
    def ref_species = (meta.reference_species_id instanceof List)
        ? meta.reference_species_id.join('').trim()
        : (meta.reference_species_id ?: '').toString().trim()
    // Defensive cleanup for samplesheets that already baked in a malformed name
    // (e.g. 'Blachea .' from a 'spp.' nominal id): collapse whitespace and drop a
    // trailing punctuation-only token, which NCBI rejects ('No such species').
    def query_species = (ref_species ?: meta.nominal_species_id).toString().replaceAll(/\s+/, ' ').replaceAll(/\s+[^A-Za-z]+$/, '').trim()
    def effective_args = ["--species '${query_species}'", "--outfolder .", args].findAll { it?.trim() }.join(' ')
    """
    # Per-task matplotlib cache + headless backend so findMitoReference does not
    # warn/fail on a read-only HOME inside the container.
    export MPLBACKEND=Agg
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    mkdir -p .mplconfig
    # Polite default delay between NCBI calls (maxForks=1 already serialises them).
    export NCBI_DELAY="\${NCBI_DELAY:-0.35}"

    findMitoReference.py \\
        ${ncbi_api_key} \\
        --species "${query_species}" \\
        --outfolder . \\
        $args

    cat <<-END_TOOL_PARAMS > 04_mitohifi_findmitoreference.tool_params_mqcrow.html
    <tr><td>MitoHiFi Find Reference</td><td><samp>${effective_args}</samp></td><td>Finds a mitochondrial reference for ${meta.id} using query species ${query_species} (nominal: ${meta.nominal_species_id}).</td></tr>
    END_TOOL_PARAMS

    # The mitohifi CLI mis-reports its version, so pin to the container tag.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: 3.2.3
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    // Older samplesheets lack the reference_species_id column; nf-schema fills the
    // missing meta value with an empty list (not null), so coerce to a clean string
    // before testing it and fall back to the raw nominal name when it's empty.
    def ref_species = (meta.reference_species_id instanceof List)
        ? meta.reference_species_id.join('').trim()
        : (meta.reference_species_id ?: '').toString().trim()
    // Defensive cleanup for samplesheets that already baked in a malformed name
    // (e.g. 'Blachea .' from a 'spp.' nominal id): collapse whitespace and drop a
    // trailing punctuation-only token, which NCBI rejects ('No such species').
    def query_species = (ref_species ?: meta.nominal_species_id).toString().replaceAll(/\s+/, ' ').replaceAll(/\s+[^A-Za-z]+$/, '').trim()
    def effective_args = ["--species '${query_species}'", "--outfolder .", args].findAll { it?.trim() }.join(' ')
    """
    touch accession.fasta
    touch accession.gb

    cat <<-END_TOOL_PARAMS > 04_mitohifi_findmitoreference.tool_params_mqcrow.html
    <tr><td>MitoHiFi Find Reference</td><td><samp>${effective_args}</samp></td><td>Finds a mitochondrial reference for ${meta.id} using query species ${query_species} (nominal: ${meta.nominal_species_id}).</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: 3.2.3
    END_VERSIONS
    """
}
