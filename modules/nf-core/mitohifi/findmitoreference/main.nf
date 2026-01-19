process MITOHIFI_FINDMITOREFERENCE {
    tag "${meta.id} $species"
    label 'process_single'
    label 'error_retry'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:master'

    input:
    tuple val(meta)

    output:
    tuple val(meta), path("*.fasta"), path("*.gb")  , emit: reference
    tuple val(meta), path("versions.yml")           , emit: versions_tuple
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    findMitoReference.py \\
        --species "$meta.nominal_species_id" \\
        --outfolder . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
    """

    stub:
    """
    touch accession.fasta
    touch accession.gb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
    END_VERSIONS
    """
}
