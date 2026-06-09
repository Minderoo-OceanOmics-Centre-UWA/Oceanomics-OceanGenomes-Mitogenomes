process MITOHIFI_FINDMITOREFERENCE {
    tag "${meta.id} $species"
    label 'process_single'
    label 'error_retry'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:master'

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
    def args = task.ext.args ?: ''
    // Prefer the NCBI-valid name resolved at samplesheet creation; fall back to
    // the raw nominal name for older samplesheets that lack the column. The raw
    // nominal name is frequently rejected by NCBI ('No such species in NCBI!').
    def query_species = (meta.reference_species_id?.trim()) ? meta.reference_species_id : meta.nominal_species_id
    def effective_args = ["--species '${query_species}'", "--outfolder .", args].findAll { it?.trim() }.join(' ')
    """
    findMitoReference.py \\
        --species "${query_species}" \\
        --outfolder . \\
        $args

    cat <<-END_TOOL_PARAMS > 04_mitohifi_findmitoreference.tool_params_mqcrow.html
    <tr><td>MitoHiFi Find Reference</td><td><samp>${effective_args}</samp></td><td>Finds a mitochondrial reference for ${meta.id} using query species ${query_species} (nominal: ${meta.nominal_species_id}).</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def query_species = (meta.reference_species_id?.trim()) ? meta.reference_species_id : meta.nominal_species_id
    def effective_args = ["--species '${query_species}'", "--outfolder .", args].findAll { it?.trim() }.join(' ')
    """
    touch accession.fasta
    touch accession.gb

    cat <<-END_TOOL_PARAMS > 04_mitohifi_findmitoreference.tool_params_mqcrow.html
    <tr><td>MitoHiFi Find Reference</td><td><samp>${effective_args}</samp></td><td>Finds a mitochondrial reference for ${meta.id} using query species ${query_species} (nominal: ${meta.nominal_species_id}).</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
    """
}
