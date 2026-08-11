process ALLOCATE_ENA_LOCUS_TAGS {
    tag "$meta.full_seqid"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(sample_fa), path(sample_tbl)
    path db_config

    output:
    tuple val(meta), path("tagged/*.tbl"), emit: tagged_tbl
    tuple val(meta), path("*.locus_tag_mapping.tsv"), emit: mapping
    tuple val(meta), path("18_ena_locus_tags.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Resolved from the candidate's technology in the subworkflow. Missing means
    // the meta never passed through EnaTargets, so refuse rather than guess a
    // prefix and write tags belonging to another study's namespace.
    def prefix = meta.ena_locus_prefix?.toString()?.trim()
    if (!prefix) {
        error "meta.ena_locus_prefix is not set for ${meta.full_seqid}: cannot allocate ENA locus tags"
    }
    """
    mkdir -p tagged
    allocate_ena_locus_tags.py \
        --og-id '${meta.id}' \
        --full-seqid '${meta.full_seqid}' \
        --input-tbl '${sample_tbl}' \
        --output-tbl 'tagged/${meta.full_seqid}.tbl' \
        --mapping '${meta.full_seqid}.locus_tag_mapping.tsv' \
        --locus-prefix '${prefix}' \
        --db-config '${db_config}'

    printf '%s\n' '<tr><td>ENA locus tags</td><td><samp>${prefix}_${String.format("%06d", meta.id.substring(2) as int)}NNN</samp></td><td>Allocates stable specimen-derived locus tags and injects them into the feature table for ${meta.full_seqid}.</td></tr>' > 18_ena_locus_tags.tool_params_mqcrow.html
    printf '"%s":\n    python: "%s"\n    allocate_ena_locus_tags: "1.0.0"\n' \
        "${task.process}" "\$(python --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    mkdir -p tagged
    cp '${sample_tbl}' 'tagged/${meta.full_seqid}.tbl'
    printf 'full_seqid\tfeature_key\tfeature_type\tcanonical_gene\tgene_occurrence\tstart\tend\tstrand\tlocus_tag\n' > '${meta.full_seqid}.locus_tag_mapping.tsv'
    printf '%s\n' '<tr><td>ENA locus tags</td><td><samp>stub</samp></td><td>Stub locus-tag allocation for ${meta.full_seqid}.</td></tr>' > 18_ena_locus_tags.tool_params_mqcrow.html
    printf '"%s":\n    python: "stub"\n    allocate_ena_locus_tags: "stub"\n' "${task.process}" > versions.yml
    """
}
