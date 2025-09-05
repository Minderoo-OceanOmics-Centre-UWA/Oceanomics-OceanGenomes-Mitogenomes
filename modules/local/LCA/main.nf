process LCA {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/lca:0.2' :
        'tylerpeirce/lca:0.2' }"

    input:
    tuple val(meta), path(blast_filtered), val(gene_type), val(annotation_name)
    path(taxdump)

    output:
    tuple val(meta), path("lca.${gene_type}.${annotation_name}.tsv"), emit: lca
    path "versions_LCA.yml", emit: versions

    script:
    """          
    export TAXONKIT_DB=$taxdump

    python /computeLCA.py \\
        $blast_filtered \\
        --output lca.${gene_type}.${annotation_name}.tsv \\
        --gene-type ${gene_type}

    
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