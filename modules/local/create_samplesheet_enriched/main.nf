process CREATE_SAMPLESHEET_ENRICHED {
    tag "Creating enriched samplesheet from ${input_files}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    path(input_files)
    val output_name
    path sql_config

    output:
    path "${output_name}", emit: samplesheet
    path "versions.yml", emit: versions

    script:
    """
    create_samplesheet.py \
        --output "${output_name}" \
        --sql-config "${sql_config}" \
        --input-files ${input_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9"
        psycopg2: "2.9.5"
    END_VERSIONS
    """

    stub:
    """
    : > ${output_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        psycopg2: "stub"
    END_VERSIONS
    """
}
