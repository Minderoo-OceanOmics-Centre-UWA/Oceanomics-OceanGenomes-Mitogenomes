process PREPARE_ENA_METADATA {
    tag "$meta.full_seqid"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    val meta
    path db_config

    output:
    tuple val(meta), path("*.ena_input_metadata.json"), emit: metadata
    path "versions.yml", emit: versions

    script:
    // The manifest STUDY line is built from this metadata, so a wrong or missing
    // value here is data landing in the wrong ENA study.
    def ena_study = meta.ena_study?.toString()?.trim()
    if (!ena_study) {
        error "meta.ena_study is not set for ${meta.full_seqid}: cannot prepare ENA metadata"
    }
    """
    prepare_ena_metadata.py \
        --config '${db_config}' \
        --og-id '${meta.id}' \
        --assembly-prefix '${meta.mt_assembly_prefix}' \
        --annotation-version '${meta.annotation_version}' \
        --full-seqid '${meta.full_seqid}' \
        --tech '${meta.sequencing_type}' \
        --seq-date '${meta.date}' \
        --code '${meta.mt_assembly_prefix.tokenize(".")[3]}' \
        --study '${ena_study}' \
        --scientific-name '${meta.scientific_name}' \
        --output '${meta.full_seqid}.ena_input_metadata.json'
    printf '"%s":\n    python: "%s"\n    prepare_ena_metadata: "1.0.0"\n' \
        "${task.process}" "\$(python --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf '{"schema_version":1,"og_id":"${meta.id}","assembly_prefix":"${meta.mt_assembly_prefix}","annotation_version":"${meta.annotation_version}","full_seqid":"${meta.full_seqid}","study":"${meta.ena_study}","biosample_accession":"SAMEA1","biosample_source":"stub","mean_depth":100,"program":"stub 1.0","platform":"${meta.sequencing_type == "hifi" ? "PACBIO_SMRT" : "ILLUMINA"}","scientific_name":"${meta.scientific_name}","run_accessions":[]}\n' > '${meta.full_seqid}.ena_input_metadata.json'
    printf '"%s":\n    python: "stub"\n    prepare_ena_metadata: "stub"\n' "${task.process}" > versions.yml
    """
}
