// Write the annotation statistics into mitogenome_data's annotation columns.
//
// Pure pusher: it takes the CSV that ANNOTATION_STATS already computed and does
// nothing but the upsert. The statistics half used to live here too, which coupled
// every gene-count verdict to having a database -- see modules/local/annotation_stats.
process PUSH_MTDNA_ANNOTATION_RESULTS {
    tag "$meta.id"
    label 'process_upload'
    
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(annotation_stats_csv)
    path config

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.annotation.upload.txt"), emit: upload
    tuple val(meta), path("12_push_mtdna_annotation_results.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args2 = task.ext.args2 ?: ''
    def effective_args = "push_emma_annotation_results.py ${args2} ${config} ${meta.id} ${annotation_stats_csv}"
    """
    # Push the results to SQL database
    push_emma_annotation_results.py \\
        $args2 \\
        $config \\
        ${meta.id} \\
        ${annotation_stats_csv} \\
        > ${meta.mt_assembly_prefix}.annotation.upload.txt

    cat <<-END_TOOL_PARAMS > 12_push_mtdna_annotation_results.tool_params_mqcrow.html
    <tr><td>Push mtDNA Annotation Results</td><td><samp>${effective_args}</samp></td><td>Uploads annotation results for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        push_emma_annotation_results: "1.0.0"
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def effective_args = "push_emma_annotation_results.py ${args2} ${config} ${meta.id} ${annotation_stats_csv}"
    """
    touch ${meta.mt_assembly_prefix}.annotation.upload.txt

    cat <<-END_TOOL_PARAMS > 12_push_mtdna_annotation_results.tool_params_mqcrow.html
    <tr><td>Push mtDNA Annotation Results</td><td><samp>${effective_args}</samp></td><td>Uploads annotation results for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        push_emma_annotation_results: "1.0.0"
    END_VERSIONS
    """
}
