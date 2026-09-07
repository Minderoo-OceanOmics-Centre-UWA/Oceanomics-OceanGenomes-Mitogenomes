process PUSH_MTDNA_ANNOTATION_RESULTS {
    tag "$meta.id"
    label 'process_upload'
    
    conda "conda-forge::python=3.9 conda-forge::psycopg2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(annotations), path(lca_combined)
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
    // Pass the resolved genetic code so annotation_stats applies the right
    // completeness profile: only code 2 is judged against the vertebrate 37-gene
    // set and gene order. Everything else (cnidarians encode only ~2 mt tRNAs;
    // no invertebrate follows the vertebrate order) is judged on the conserved
    // PCG+rRNA core. The class is still passed as the fallback selector for a
    // meta with no resolved code.
    def class_arg = meta.class ? "--class '${meta.class}'" : ''
    def gcode_arg = meta.genetic_code ? "--genetic-code ${meta.genetic_code}" : ''
    // Taxonomy for the curated gene-order variant lookup ONLY. It never affects
    // completeness, only which non-canonical orders are accepted for this clade.
    // Same shape as modules/local/reference_divergence/main.nf, which is the
    // established pattern for reading taxonomy off meta.
    def family      = (meta.family ?: '').toString().trim()
    def taxon_order = (meta.order ?: '').toString().trim()
    // There is no separate genus field on meta; the first whitespace token of the
    // nominal species id is how bin/reference_divergence_check.py derives it too.
    def genus       = (meta.nominal_species_id ?: '').toString().trim().split(/\s+/)[0] ?: ''
    def family_arg  = family ? "--family '${family}'" : ''
    def order_arg   = taxon_order ? "--order '${taxon_order}'" : ''
    def genus_arg   = genus ? "--genus '${genus}'" : ''
    // The second taxonomy opinion, for the advisory order_variant_taxon_check.
    // Optional: an assembly with no lca_combined (the file is staged as a
    // placeholder) simply records 'no'.
    def lca_arg     = lca_combined ? "--lca-combined ${lca_combined}" : ''
    def effective_args = ["annotation_stats.py ${args} ${class_arg} ${gcode_arg} ${family_arg} ${order_arg} ${genus_arg} ${lca_arg} *.gff proteins", "push_emma_annotation_results.py ${args2} ${config} ${meta.id} ${meta.mt_assembly_prefix}.annotation_stats.csv"].findAll { it?.trim() }.join('; ')
    """
    # Compile the statistics.
    #
    # imports mito_gene_order.py and orf_utils.py -- named here on purpose. Nextflow
    # hashes only those bin/ scripts whose filenames appear as tokens in a task's
    # command script, and neither of these is ever invoked directly: they reach the
    # task as a Python import. Without this line, editing the shared gene-order table
    # alone and resuming would re-run NOTHING and silently return stale results.
    annotation_stats.py \\
        $args \\
        ${class_arg} \\
        ${gcode_arg} \\
        ${family_arg} \\
        ${order_arg} \\
        ${genus_arg} \\
        ${lca_arg} \\
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
