process LCA {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0' :
        'quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0' }"

    input:
    tuple val(meta), path(blast), val(gene_type), val(annotation_name)
    path(worms)

    output:
    tuple val(meta), path("lca.${gene_type}.${annotation_name}.tsv"), emit: lca
    path("lca_raw.${gene_type}.${annotation_name}.tsv"), emit: lca_raw
    path("lca_short.${gene_type}.${annotation_name}.tsv"), emit: lca_short
    tuple val(meta), path("09_lca.${gene_type}.${annotation_name}.tool_params_mqcrow.html"), emit: tool_params
    path "versions_LCA.yml", emit: versions

    script:
    def effective_args = "--file ${blast} --output lca_short.${gene_type}.${annotation_name}.tsv --worms_file ${worms} --raw_output lca_raw.${gene_type}.${annotation_name}.tsv --final_output lca.${gene_type}.${annotation_name}.tsv --seq_type ${gene_type}"
    """          
    calculateLCA.py \\
        --file $blast \\
        --output lca_short.${gene_type}.${annotation_name}.tsv \\
        --worms_file $worms \\
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
    def effective_args = "--file ${blast} --output lca_short.${gene_type}.${annotation_name}.tsv --worms_file ${worms} --raw_output lca_raw.${gene_type}.${annotation_name}.tsv --final_output lca.${gene_type}.${annotation_name}.tsv --seq_type ${gene_type}"
    """
    touch lca.${gene_type}.${annotation_name}.tsv

    cat <<-END_TOOL_PARAMS > 09_lca.${gene_type}.${annotation_name}.tool_params_mqcrow.html
    <tr><td>LCA</td><td><samp>${effective_args}</samp></td><td>Calculates lowest common ancestor assignments for ${gene_type} sequence ${annotation_name} from filtered BLAST results.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions_LCA.yml
    "${task.process}":
        Python: \$(python -V | sed 's/Python //g')
    END_VERSIONS
    """
}
