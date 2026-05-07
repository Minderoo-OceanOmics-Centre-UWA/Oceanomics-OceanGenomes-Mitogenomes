process SPECIES_VALIDATION {
    tag "$meta.id"
    label 'process_medium'
    
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(blast_files), path(lca_files)
    path config

    output:
    tuple val(meta), path("lca_results.${meta.id}.tsv"), emit: summary
    tuple val(meta), path("lca_combined.${meta.id}.tsv"), path("blast_combined.${meta.id}.tsv"), emit: full
    tuple val(meta), path("11_species_validation.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def lca_files_str = lca_files instanceof List ? lca_files.join(',') : lca_files
    def blast_files_str = blast_files instanceof List ? blast_files.join(',') : blast_files
    def effective_args = [args, config, meta.id, lca_files_str, blast_files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    species_validation.py \\
        $args \\
        $config \\
        ${meta.id} \\
        "${lca_files_str}" \\
        "${blast_files_str}"

    cat <<-END_TOOL_PARAMS > 11_species_validation.tool_params_mqcrow.html
    <tr><td>Species Validation</td><td><samp>${effective_args}</samp></td><td>Compares filtered BLAST/LCA calls against the nominal species for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        species_validation: "1.0.0"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def lca_files_str = lca_files instanceof List ? lca_files.join(',') : lca_files
    def blast_files_str = blast_files instanceof List ? blast_files.join(',') : blast_files
    def effective_args = [args, config, meta.id, lca_files_str, blast_files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    touch lca_results.${meta.id}.tsv
    touch lca_combined.${meta.id}.tsv
    touch blast_combined.${meta.id}.tsv

    cat <<-END_TOOL_PARAMS > 11_species_validation.tool_params_mqcrow.html
    <tr><td>Species Validation</td><td><samp>${effective_args}</samp></td><td>Compares filtered BLAST/LCA calls against the nominal species for ${meta.id}.</td></tr>
    END_TOOL_PARAMS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        species_validation: "1.0.0"
    END_VERSIONS
    """
}
