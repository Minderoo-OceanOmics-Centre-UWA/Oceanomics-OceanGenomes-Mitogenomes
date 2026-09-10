// Write the lca_validation row from the record SPECIES_VALIDATION produced.
//
// Pure pusher. The species comparison that decides what goes in the row is
// DB-free and runs unconditionally (see modules/local/species_validation); this
// process is the part that needs a database, and is therefore the part that is
// skipped when uploads are off.
process PUSH_SPECIES_VALIDATION {
    tag "$meta.id"
    label 'process_upload'

    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(validation_record)
    path config

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix ?: meta.id}.species_validation.upload.txt"), emit: upload
    tuple val(meta), path("11_push_species_validation.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    def effective_args = [args, config, validation_record].findAll { it?.toString()?.trim() }.join(' ')
    """
    push_species_validation.py \\
        $args \\
        $config \\
        $validation_record \\
        > ${mt_assembly_prefix}.species_validation.upload.txt

    cat <<-END_TOOL_PARAMS > 11_push_species_validation.tool_params_mqcrow.html
    <tr><td>Push Species Validation</td><td><samp>push_species_validation.py ${effective_args}</samp></td><td>Upserts the lca_validation row for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        push_species_validation: "1.0.0"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    def effective_args = [args, config, validation_record].findAll { it?.toString()?.trim() }.join(' ')
    """
    touch ${mt_assembly_prefix}.species_validation.upload.txt

    cat <<-END_TOOL_PARAMS > 11_push_species_validation.tool_params_mqcrow.html
    <tr><td>Push Species Validation</td><td><samp>push_species_validation.py ${effective_args}</samp></td><td>Upserts the lca_validation row for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        push_species_validation: "1.0.0"
    END_VERSIONS
    """
}
