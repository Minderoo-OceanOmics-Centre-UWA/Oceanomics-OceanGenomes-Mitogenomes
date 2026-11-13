process BUILD_SOURCE_MODIFIERS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    val meta
    path config

    output:
    tuple val(meta), path("${meta.id}.bankit_metadata.csv")                , emit: bankit_metadata
    tuple val(meta), path("${meta.id}.bankit_metadata_latlon_cleaned.csv") , emit: cleaned_metadata
    // Named exactly, and NOT optional.
    //
    // It was `path("${meta.mt_assembly_prefix}.*.src"), optional: true`, which let the script
    // write .src files under OTHER variants' names, match none of them, and exit 0. The sample
    // was then dropped by the inner join in ena_submission_prep and reported as a held row
    // blaming a task failure that never happened. An exact required name turns that into a
    // task failure that says which file is missing, which the global errorStrategy 'ignore'
    // still confines to this one assembly.
    tuple val(meta), path("${meta.full_seqid}.src")                        , emit: src_file
    tuple val(meta), path("16_build_source_modifiers.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def effective_args = "--config ${config} --og-id ${meta.id} --seq-tech ${meta.sequencing_type} --assembly-id ${meta.mt_assembly_prefix} --full-seqid ${meta.full_seqid}"
    """

    build_source_modifiers.py \\
        --config $config \\
        --og-id "${meta.id}" \\
        --seq-tech "${meta.sequencing_type}" \\
        --assembly-id "${meta.mt_assembly_prefix}" \\
        --full-seqid "${meta.full_seqid}"

    cat <<-END_TOOL_PARAMS > 16_build_source_modifiers.tool_params_mqcrow.html
    <tr><td>Build Source Modifiers</td><td><samp>${effective_args}</samp></td><td>Builds GenBank source modifier files for ${meta.id} from SQL metadata.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def effective_args = "--config ${config} --og-id ${meta.id} --seq-tech ${meta.sequencing_type} --assembly-id ${meta.mt_assembly_prefix} --full-seqid ${meta.full_seqid}"
    """
    # Names must match the output: block exactly. src_file is no longer optional, so a
    # mismatch here fails the stub loudly instead of silently truncating every -stub run of
    # ENA_SUBMISSION_PREP before table2asn.
    : > ${meta.id}.bankit_metadata.csv
    : > ${meta.id}.bankit_metadata_latlon_cleaned.csv
    : > ${meta.full_seqid}.src

    cat <<-END_TOOL_PARAMS > 16_build_source_modifiers.tool_params_mqcrow.html
    <tr><td>Build Source Modifiers</td><td><samp>${effective_args}</samp></td><td>Builds GenBank source modifier files for ${meta.id} from SQL metadata.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //' 2>/dev/null || echo "3.9.0")
    END_VERSIONS
    """
}
