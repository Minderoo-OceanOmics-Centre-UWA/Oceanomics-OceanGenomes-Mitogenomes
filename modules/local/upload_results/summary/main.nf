process UPLOAD_RESULTS_SUMMARY {
    tag "upload_results_summary"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path upload_files, stageAs: "upload_inputs/*"

    output:
    path "upload_results_summary_mqc.tsv", emit: summary
    path "upload_results_appendix.txt",    emit: appendix
    path "14_upload_results_summary.tool_params_mqcrow.html", emit: tool_params
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def effective_args = "--input-dir upload_inputs --summary upload_results_summary_mqc.tsv --appendix upload_results_appendix.txt"
    """
    compile_upload_report.py \\
        --input-dir upload_inputs \\
        --summary upload_results_summary_mqc.tsv \\
        --appendix upload_results_appendix.txt

    cat <<-END_TOOL_PARAMS > 14_upload_results_summary.tool_params_mqcrow.html
    <tr><td>Upload Results Summary</td><td><samp>${effective_args}</samp></td><td>Consolidates per-step SQL upload status files into a single summary TSV and appendix.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        compile_upload_report: "1.0.0"
    END_VERSIONS
    """

    stub:
    def effective_args = "--input-dir upload_inputs --summary upload_results_summary_mqc.tsv --appendix upload_results_appendix.txt"
    """
    touch upload_results_summary_mqc.tsv
    touch upload_results_appendix.txt

    cat <<-END_TOOL_PARAMS > 14_upload_results_summary.tool_params_mqcrow.html
    <tr><td>Upload Results Summary</td><td><samp>${effective_args}</samp></td><td>Consolidates per-step SQL upload status files into a single summary TSV and appendix.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        compile_upload_report: "stub"
    END_VERSIONS
    """
}
