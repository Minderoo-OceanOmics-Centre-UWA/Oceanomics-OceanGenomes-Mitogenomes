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
    // Non-validating diagnostic, published alongside the authoritative outputs.
    // Nothing consumes it as a channel: PUSH_SPECIES_VALIDATION writes
    // lca_validation from validation_record above, and this has no path there.
    tuple val(meta), path("taxonomic_evidence.${meta.mt_assembly_prefix ?: meta.id}.tsv"), emit: taxonomic_evidence
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
    // Resolved lineage from the samplesheet, so a label above family can be pinned to
    // the rank it actually claims. No Latin suffix separates class from order
    // (Asteroidea/Scleralcyonacea are both -acea/-oidea), so without this
    // 'Asteroidea sp.' reads as a genus that does not exist. Reporting only: it never
    // changes Found_in_blast_YN. Passed only when present, like nominal_arg.
    def lineage_pairs = [class: meta.class, order: meta.order, family: meta.family]
        .findAll { _k, v -> v?.toString()?.trim() }
        .collect { k, v -> "${k}=${v.toString().trim()}" }
        .join(',')
    def lineage_arg = lineage_pairs ? "--lineage '${lineage_pairs}'" : ''
    def effective_args = [args, nominal_arg, lineage_arg, "--assembly-prefix ${mt_assembly_prefix}", meta.id, lca_files_str, blast_files_str].findAll { it?.toString()?.trim() }.join(' ')
    """
    species_validation.py \\
        $args \\
        ${nominal_arg} \\
        ${lineage_arg} \\
        --assembly-prefix ${mt_assembly_prefix} \\
        --record-file validation_record.${mt_assembly_prefix}.json \\
        ${meta.id} \\
        "${lca_files_str}" \\
        "${blast_files_str}"

    # Non-validating evidence view. Derived ONLY from the files species_validation.py
    # just wrote; it has no database connection and no write path to
    # lca_validation, Found_in_blast_YN, proceed_qc or submission metadata. It
    # exists because "not validated (no Found_in_blast_YN=Yes)" was the identical
    # reason given for 68/68 assemblies in NOVA_260909_LA, covering three quite
    # different situations that need three different actions.
    #
    # --validated is read back OUT of the authoritative result so the report can
    # agree with it; this step never sets it.
    # Read the authoritative verdict back out of the record species_validation.py
    # just wrote: action=upsert is exactly "validated, ask for an lca_validation
    # row". Do NOT scrape lca_results.tsv for a trailing Yes -- validated_rank is
    # appended after Found_in_blast_YN, so Yes is no longer the last column.
    validated_flag=""
    if grep -q '"action": *"upsert"' validation_record.${mt_assembly_prefix}.json 2>/dev/null; then
        validated_flag="--validated"
    fi

    taxonomic_evidence.py \\
        --assembly "${mt_assembly_prefix}" \\
        --og-id "${meta.id}" \\
        --nominal "${meta.nominal_species_id ?: ''}" \\
        --blast "blast_combined.${mt_assembly_prefix}.tsv" \\
        --lca "lca_combined.${mt_assembly_prefix}.tsv" \\
        --class "${meta.class ?: ''}" \\
        --family "${meta.family ?: ''}" \\
        --order "${meta.order ?: ''}" \\
        \$validated_flag \\
        --out "taxonomic_evidence.${mt_assembly_prefix}.tsv"

    cat <<-END_TOOL_PARAMS > 11_species_validation.tool_params_mqcrow.html
    <tr><td>Species Validation</td><td><samp>${effective_args}</samp></td><td>Compares filtered BLAST/LCA calls against the nominal species for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        species_validation: "1.0.0"
        taxonomic_evidence: "1.0.0"
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
    touch taxonomic_evidence.${mt_assembly_prefix}.tsv

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
