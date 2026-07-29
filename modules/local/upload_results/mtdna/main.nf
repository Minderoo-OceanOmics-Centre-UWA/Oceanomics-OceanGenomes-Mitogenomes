process PUSH_MTDNA_ASSM_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    // depth_tsv is the uniform remap-based depth from MITOGENOME_COVERAGE, or the
    // header-only assets/empty_mito_depth.tsv placeholder for a row that never
    // reached annotation (failed, under-length, or a discarded assembly variant).
    tuple val(meta), path(fasta), path(out_log), path(depth_tsv)
    path config

    output:
    path "${meta.mt_assembly_prefix ?: meta.id}.mtdna.upload.txt", emit: upload
    tuple val(meta), path("10_push_mtdna_assm_results.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    // Corrected circular verdict from the assembly subworkflow (GetOrganelle
    // reference check / MitoHiFi circularity check); null/absent -> no override.
    def circular = meta.circular == null ? 'null' : meta.circular.toString()
    def depth_arg = depth_tsv ? "--depth-tsv ${depth_tsv}" : ''
    def effective_args = [args, "--circular ${circular}", depth_arg, config, mt_assembly_prefix, out_log, fasta].findAll { it?.toString()?.trim() }.join(' ')
    """
    push_mtdna_assm_results.py \\
        $args \\
        --circular ${circular} \\
        $depth_arg \\
        $config \\
        ${mt_assembly_prefix} \\
        $out_log \\
        $fasta \\
        > ${mt_assembly_prefix}.mtdna.upload.txt

    cat <<-END_TOOL_PARAMS > 10_push_mtdna_assm_results.tool_params_mqcrow.html
    <tr><td>Push mtDNA Assembly Results</td><td><samp>${effective_args}</samp></td><td>Uploads mitogenome assembly statistics for ${meta.id} to the SQL database.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        push_mtdna_assm_results: "1.0.0"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    def effective_args = [args, config, mt_assembly_prefix, out_log, fasta, depth_tsv].findAll { it?.toString()?.trim() }.join(' ')
    """
    touch ${mt_assembly_prefix}.mtdna.upload.txt

    cat <<-END_TOOL_PARAMS > 10_push_mtdna_assm_results.tool_params_mqcrow.html
    <tr><td>Push mtDNA Assembly Results</td><td><samp>${effective_args}</samp></td><td>Uploads mitogenome assembly statistics for ${meta.id} to the SQL database.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        push_mtdna_assm_results: "1.0.0"
    END_VERSIONS
    """
}
