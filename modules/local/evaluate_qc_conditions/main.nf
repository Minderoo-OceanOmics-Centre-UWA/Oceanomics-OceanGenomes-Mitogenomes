process EVALUATE_QC_CONDITIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(blast_table), path(annotation_csv), path(circularity_check)

    output:
    tuple val(meta), path("species_name.txt"), path("proceed_qc.txt"), path("circular.txt"), emit: evaluation
    tuple val(meta), path("14_evaluate_qc_conditions.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // meta.circular is the GetOrganelle topology verdict (true/false/null). MitoHiFi
    // leaves it null; for those samples the verdict comes from the circularity-check
    // sidecar (final_verdict_circular) inside the script.
    def circular_flag = (meta.circular == null) ? 'null' : meta.circular.toString()
    def effective_args = "--blast-table ${blast_table} --annotation-csv ${annotation_csv} --circularity-check ${circularity_check} --circular ${circular_flag} --output-species species_name.txt --output-proceed proceed_qc.txt --output-circular circular.txt --output-versions versions.yml"
    """
    evaluate_qc_conditions.py \\
        --blast-table ${blast_table} \\
        --annotation-csv ${annotation_csv} \\
        --circularity-check ${circularity_check} \\
        --circular ${circular_flag} \\
        --output-species species_name.txt \\
        --output-proceed proceed_qc.txt \\
        --output-circular circular.txt \\
        --output-versions versions.yml

    cat <<-END_TOOL_PARAMS > 14_evaluate_qc_conditions.tool_params_mqcrow.html
    <tr><td>Evaluate QC Conditions</td><td><samp>${effective_args}</samp></td><td>Checks species validation and annotation QC conditions for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    """

    stub:
    def effective_args = "--blast-table ${blast_table} --annotation-csv ${annotation_csv} --output-species species_name.txt --output-proceed proceed_qc.txt --output-circular circular.txt --output-versions versions.yml"
    """
    echo "test_species" > species_name.txt
    echo "true" > proceed_qc.txt
    echo "true" > circular.txt

    cat <<-END_TOOL_PARAMS > 14_evaluate_qc_conditions.tool_params_mqcrow.html
    <tr><td>Evaluate QC Conditions</td><td><samp>${effective_args}</samp></td><td>Checks species validation and annotation QC conditions for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
