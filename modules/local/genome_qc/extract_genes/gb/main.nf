process EXTRACT_GENES_GB {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(tbl)


    output:
    tuple val(meta), path('cds'), emit: cds_dir
    tuple val(meta), path('cds/*.fa'), emit: cds_fastas
    tuple val(meta), path('genes'), optional: true, emit: genes_dir
    tuple val(meta), path('genes/*.fa'), optional: true, emit: genes_fastas
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when ?: true


    script:
    def asm = (meta.mt_assembly_prefix ?: meta.sample_id ?: fasta.baseName)
    // Per-sample mitochondrial translation table, same idiom as GEN_FILES_TABLE2ASN
    // and FORMAT_FILES. It becomes the [mgcode=] tag on every extracted CDS header.
    def gcode = task.ext.code ?: meta.genetic_code

    """
    extract_cds_from_tbl.py \\
        --fasta ${fasta} \\
        --tbl ${tbl} \\
        --outdir . \\
        --assembly ${asm} \\
        --genetic-code ${gcode}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    asm=${meta.mt_assembly_prefix ?: (meta.sample_id ?: "stub")}
    mkdir -p cds genes
    : > cds/\${asm}.fa
    # genes output is optional; create directory to satisfy path output
    : > genes/\${asm}.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
