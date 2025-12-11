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
    path "versions_LCA.yml", emit: versions

    script:
    """          
    calculateLCA.py \\
        --file $blast \\
        --output lca_short.${gene_type}.${annotation_name}.tsv \\
        --worms_file $worms \\
        --raw_output lca_raw.${gene_type}.${annotation_name}.tsv \\
        --final_output lca.${gene_type}.${annotation_name}.tsv \\
        --seq_type ${gene_type}
     
    cat <<-END_VERSIONS > versions_LCA.yml
    "${task.process}":
        Python: \$(python -V | sed 's/Python //g')
        TaxonKit: \$(/taxonkit version | sed 's/taxonkit //g')
    END_VERSIONS
    
    """

    stub:
    """
    touch lca.${gene_type}.${annotation_name}.tsv

    cat <<-END_VERSIONS > versions_LCA.yml
    "${task.process}":
        Python: \$(python -V | sed 's/Python //g')
        TaxonKit: \$(/taxonkit version | sed 's/taxonkit //g')
    END_VERSIONS
    """
}