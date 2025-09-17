process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5222a42b366a0468a4c795f5057c2b8cfe39489548f8bd807e8ac0f80069bad5/data':
        'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd' }"

    input:
    tuple val(meta), path(fasta), val(gene_type), val(annotation_name)
    val(curated_blast_db)
    path(taxdb)

    output:
    tuple val(meta), path("blast.${gene_type}.${annotation_name}.filtered.tsv"), val(gene_type), val(annotation_name), emit: filtered
    tuple val(meta), path("blast.${gene_type}.${annotation_name}.filtered.tsv"), emit: validation
    path "filtered_summary.${gene_type}.txt"   , emit: summary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "blast.${gene_type}.${annotation_name}"

    """
    mkdir -p lca
    
    blastn \\
        -num_threads ${task.cpus} \\
        -db ${curated_blast_db} \\
        -query ${fasta} \\
        ${args} \\
        -out ${prefix}.tsv

    # Filter the results
    awk -F '\t' '{if ((\$15 - \$14 > 200) && (\$7 > 98)) print}' ${prefix}.tsv > ${prefix}.filtered.tsv
    
    # Add in a column for the days date and the gene type
    sed -i "s/\$/\\t\$(date +%y%m%d)/" ${prefix}.filtered.tsv
    wait
    sed -i "s/\$/\\t${gene_type}/" ${prefix}.filtered.tsv

    # Create a summary for multiqc
    awk -F '\t' '{combo[\$1"\t"\$4]++} END {print "Sample\tSpecies\tHits"; for (c in combo) print c, combo[c]}' ${prefix}.filtered.tsv > filtered_summary.${gene_type}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
