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
    tuple val(meta), path("${meta.mt_assembly_prefix ?: meta.id}.findmitoreference_status.tsv"), emit: status
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

    set +e
    findMitoReference.py \\
        ${ncbi_api_key} \\
        --species "${query_species}" \\
        --outfolder . \\
        $args
    reference_exit=\$?
    set -e

    ref_fasta=\$(find . -maxdepth 1 -name '*.fasta' -type f -size +0c | head -n1)
    ref_gb=\$(find . -maxdepth 1 -name '*.gb' -type f -size +0c | head -n1)
    if [ "\$reference_exit" -eq 0 ] && [ -n "\$ref_fasta" ] && [ -n "\$ref_gb" ]; then
        printf 'sample\tstatus\texit_code\tquery_species\n%s\tfound\t0\t%s\n' '${meta.id}' '${query_species}' > ${meta.mt_assembly_prefix ?: meta.id}.findmitoreference_status.tsv
    elif [ "\$reference_exit" -ne 0 ] && [ '${task.attempt}' -le 3 ]; then
        # Preserve technical failures for the first three attempts so the process
        # error strategy can retry them. On the fourth attempt the sample is kept
        # in the cohort with a structured lookup_error outcome below.
        exit "\$reference_exit"
    else
        # Exit 0 without a reference pair is a biological no-reference outcome.
        # A non-zero fourth attempt is an exhausted technical lookup error. Both
        # keep the sample routable, but remain distinguishable in audit output.
        rm -f -- ./*.fasta ./*.gb
        : > no_reference.fasta
        : > no_reference.gb
        if [ "\$reference_exit" -eq 0 ]; then
            reference_status=no_reference
        else
            reference_status=lookup_error
        fi
        printf 'sample\tstatus\texit_code\tquery_species\n%s\t%s\t%s\t%s\n' '${meta.id}' "\$reference_status" "\$reference_exit" '${query_species}' > ${meta.mt_assembly_prefix ?: meta.id}.findmitoreference_status.tsv
    fi

    printf '%s\n' '<tr><td>MitoHiFi Find Reference</td><td><samp>${effective_args}</samp></td><td>Finds a mitochondrial reference for ${meta.id} using query species ${query_species} (nominal: ${meta.nominal_species_id}).</td></tr>' > 04_mitohifi_findmitoreference.tool_params_mqcrow.html

    # The mitohifi CLI mis-reports its version, so pin to the container tag.
    printf '"%s":\n    mitohifi: 3.2.3\n' '${task.process}' > versions.yml
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
    printf 'sample\tstatus\texit_code\tquery_species\n${meta.id}\tfound\t0\t${query_species}\n' > ${meta.mt_assembly_prefix ?: meta.id}.findmitoreference_status.tsv

    printf '%s\n' '<tr><td>MitoHiFi Find Reference</td><td><samp>${effective_args}</samp></td><td>Finds a mitochondrial reference for ${meta.id} using query species ${query_species} (nominal: ${meta.nominal_species_id}).</td></tr>' > 04_mitohifi_findmitoreference.tool_params_mqcrow.html
    printf '"%s":\n    mitohifi: 3.2.3\n' '${task.process}' > versions.yml
    """
}
