// Compare the per-region LCA/BLAST calls against the sample's nominal species.
//
// Pure QC and deliberately DB-free. It produces Found_in_blast_YN, which is one of
// the two conditions the whole QC gate turns on, so binding it to a database meant
// no verdict could be reached with --skip_upload_results or without --sql_config.
// The nominal species now comes from the samplesheet via meta.nominal_species_id,
// and the lca_validation row is written by PUSH_SPECIES_VALIDATION from the record
// emitted here.
process SPECIES_VALIDATION {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(blast_files), path(lca_files)

    output:
    // Qualify these with mt_assembly_prefix (not just meta.id): an OG can have
    // multiple assembly attempts (e.g. hifi + hic + getorg reseed), and id-only
    // names collide (silent last-write-wins overwrite) when collected flat into
    // the shared species_validation/ output directory.
    tuple val(meta), path("lca_results.${meta.mt_assembly_prefix ?: meta.id}.tsv"), emit: summary
    tuple val(meta), path("lca_combined.${meta.mt_assembly_prefix ?: meta.id}.tsv"), path("blast_combined.${meta.mt_assembly_prefix ?: meta.id}.tsv"), emit: full
    tuple val(meta), path("validation_record.${meta.mt_assembly_prefix ?: meta.id}.json"), emit: validation_record
    tuple val(meta), path("11_species_validation.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    def lca_files_str = lca_files instanceof List ? lca_files.join(',') : lca_files
    def blast_files_str = blast_files instanceof List ? blast_files.join(',') : blast_files
    // The samplesheet's nominal_species_id, not a database lookup. An empty value is
    // a legitimate state (the script records the species-match columns as N/A), so
    // this is passed only when present rather than defaulted to something.
    def nominal_species = (meta.nominal_species_id ?: '').toString().trim()
    def nominal_arg = nominal_species ? "--nominal-species '${nominal_species}'" : ''
    def effective_args = [args, nominal_arg, "--assembly-prefix ${mt_assembly_prefix}", meta.id, lca_files_str, blast_files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    species_validation.py \\
        $args \\
        ${nominal_arg} \\
        --assembly-prefix ${mt_assembly_prefix} \\
        --record-file validation_record.${mt_assembly_prefix}.json \\
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
    def mt_assembly_prefix = meta.mt_assembly_prefix ?: meta.id
    def lca_files_str = lca_files instanceof List ? lca_files.join(',') : lca_files
    def blast_files_str = blast_files instanceof List ? blast_files.join(',') : blast_files
    def nominal_species = (meta.nominal_species_id ?: '').toString().trim()
    def nominal_arg = nominal_species ? "--nominal-species '${nominal_species}'" : ''
    def effective_args = [args, nominal_arg, "--assembly-prefix ${mt_assembly_prefix}", meta.id, lca_files_str, blast_files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    touch lca_results.${mt_assembly_prefix}.tsv
    touch lca_combined.${mt_assembly_prefix}.tsv
    touch blast_combined.${mt_assembly_prefix}.tsv
    printf '{"action": "skip", "reason": "stub"}\\n' > validation_record.${mt_assembly_prefix}.json

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
