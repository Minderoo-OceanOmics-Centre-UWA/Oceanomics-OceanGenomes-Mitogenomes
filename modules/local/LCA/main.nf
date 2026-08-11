process LCA {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0' :
        'quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0' }"

    input:
    tuple val(meta), path(blast), val(gene_type), val(annotation_name)
    path(worms)
    val(lca_cache_dir)
    val(lca_cache_signature)
    path(lca_script)

    output:
    // A region with no valid BLAST hits is an expected outcome for some
    // 12S/16S/CO1 queries; calculateLCA.py writes header-only lca/lca_raw files
    // in that case. These outputs are therefore NOT optional: downstream
    // grouping in upload_results_mito closes each sample's group once its
    // region count is reached, so every annotated region must yield exactly one
    // lca and one lca_raw file. Making these optional again would leave the
    // group one short and stall the sample until the end of the run.
    tuple val(meta), path("lca.${gene_type}.${annotation_name}.tsv"), emit: lca
    tuple val(meta), path("lca_raw.${gene_type}.${annotation_name}.tsv"), emit: lca_raw
    // lca_short stays optional: it is a reporting convenience, not part of the
    // per-region accounting above.
    path("lca_short.${gene_type}.${annotation_name}.tsv"), emit: lca_short, optional: true
    tuple val(meta), path("09_lca.${gene_type}.${annotation_name}.tool_params_mqcrow.html"), emit: tool_params
    path "versions_LCA.yml", emit: versions

    script:
    def effective_args = "--file ${blast} --output lca_short.${gene_type}.${annotation_name}.tsv --worms_file ${worms} --cache_dir ${lca_cache_dir} --raw_output lca_raw.${gene_type}.${annotation_name}.tsv --final_output lca.${gene_type}.${annotation_name}.tsv --seq_type ${gene_type}"
    """          
    python ${lca_script} \\
        --file $blast \\
        --output lca_short.${gene_type}.${annotation_name}.tsv \\
        --worms_file $worms \\
        --cache_dir '${lca_cache_dir}' \\
        --raw_output lca_raw.${gene_type}.${annotation_name}.tsv \\
        --final_output lca.${gene_type}.${annotation_name}.tsv \\
        --seq_type ${gene_type}

    cat <<-END_TOOL_PARAMS > 09_lca.${gene_type}.${annotation_name}.tool_params_mqcrow.html
    <tr><td>LCA</td><td><samp>${effective_args}</samp></td><td>Calculates lowest common ancestor assignments for ${gene_type} sequence ${annotation_name} from filtered BLAST results.</td></tr>
    END_TOOL_PARAMS
     
    cat <<-END_VERSIONS > versions_LCA.yml
    "${task.process}":
        Python: \$(python -V | sed 's/Python //g')
    END_VERSIONS
    
    """

    stub:
    def effective_args = "--file ${blast} --output lca_short.${gene_type}.${annotation_name}.tsv --worms_file ${worms} --cache_dir ${lca_cache_dir} --raw_output lca_raw.${gene_type}.${annotation_name}.tsv --final_output lca.${gene_type}.${annotation_name}.tsv --seq_type ${gene_type}"
    """
    touch lca.${gene_type}.${annotation_name}.tsv
    touch lca_raw.${gene_type}.${annotation_name}.tsv

    cat <<-END_TOOL_PARAMS > 09_lca.${gene_type}.${annotation_name}.tool_params_mqcrow.html
    <tr><td>LCA</td><td><samp>${effective_args}</samp></td><td>Calculates lowest common ancestor assignments for ${gene_type} sequence ${annotation_name} from filtered BLAST results.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions_LCA.yml
    "${task.process}":
        Python: \$(python -V | sed 's/Python //g')
    END_VERSIONS
    """
}
