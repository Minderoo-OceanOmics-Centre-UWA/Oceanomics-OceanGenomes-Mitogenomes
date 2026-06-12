process PUSH_MTDNA_ANNOTATION_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(annotations) 
    path config

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.annotation.upload.txt"), emit: upload
    tuple val(meta), path("${meta.mt_assembly_prefix}.annotation_stats.csv"), emit: stats
    tuple val(meta), path("12_push_mtdna_annotation_results.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // Qualify the per-assembly output names with mt_assembly_prefix (not just
    // meta.id): a sample can have multiple assemblies (e.g. hifi v321 + v323 +
    // getorg), and id-only names collide when collected flat downstream.
    def effective_args = ["annotation_stats.py ${args} *.gff proteins", "push_emma_annotation_results.py ${args2} ${config} ${meta.id} ${meta.mt_assembly_prefix}.annotation_stats.csv"].findAll { it?.trim() }.join('; ')
    """
    # Compile the statistics
    annotation_stats.py \\
        $args \\
        *.gff \\
        proteins

    wait

    # annotation_stats.py names its output from the OG id only
    # (<og_id>.annotation_stats.csv). Rename to the full assembly prefix so a
    # sample's multiple assemblies don't collide when collected flat downstream.
    mv *.annotation_stats.csv ${meta.mt_assembly_prefix}.annotation_stats.csv

    # Push the results to SQL database
    push_emma_annotation_results.py \\
        $args2 \\
        $config \\
        ${meta.id} \\
        ${meta.mt_assembly_prefix}.annotation_stats.csv \\
        > ${meta.mt_assembly_prefix}.annotation.upload.txt

    cat <<-END_TOOL_PARAMS > 12_push_mtdna_annotation_results.tool_params_mqcrow.html
    <tr><td>Push mtDNA Annotation Results</td><td><samp>${effective_args}</samp></td><td>Calculates annotation statistics and uploads annotation results for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        annotation_stats: "1.0.0"
        push_emma_annotation_results: "1.0.0"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def effective_args = ["annotation_stats.py ${args} *.gff proteins", "push_emma_annotation_results.py ${args2} ${config} ${meta.id} ${meta.mt_assembly_prefix}.annotation_stats.csv"].findAll { it?.trim() }.join('; ')
    """
    touch ${meta.mt_assembly_prefix}.annotation.upload.txt
    touch ${meta.mt_assembly_prefix}.annotation_stats.csv

    cat <<-END_TOOL_PARAMS > 12_push_mtdna_annotation_results.tool_params_mqcrow.html
    <tr><td>Push mtDNA Annotation Results</td><td><samp>${effective_args}</samp></td><td>Calculates annotation statistics and uploads annotation results for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        annotation_stats: "1.0.0"
        push_emma_annotation_results: "1.0.0"
    END_VERSIONS
    """
}
