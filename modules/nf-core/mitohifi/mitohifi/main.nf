process MITOHIFI_MITOHIFI {
    tag "$meta.id"
    label 'process_high'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    container 'ghcr.io/marcelauliano/mitohifi:master'

    input:
    tuple val(meta), path(hifi_cat), path(ref_fa), path(ref_gb)
    val input_mode
    val mito_code

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.fasta")                , emit: fasta
    tuple val(meta), path("${meta.mt_assembly_prefix}.contigs_stats.tsv")    , emit: stats
    tuple val(meta), path("${meta.mt_assembly_prefix}.gb")                   , emit: gb                         , optional: true
    // tuple val(meta), path("final_mitogenome.gff")           , emit: gff                        , optional: true
    tuple val(meta), path("all_potential_contigs.fa")       , emit: all_potential_contigs      , optional: true
    tuple val(meta), path("contigs_annotations.png")        , emit: contigs_annotations        , optional: true
    tuple val(meta), path("contigs_circularization/")       , emit: contigs_circularization    , optional: true
    tuple val(meta), path("contigs_filtering/")             , emit: contigs_filtering          , optional: true
    tuple val(meta), path("coverage_mapping/")              , emit: coverage_mapping           , optional: true
    tuple val(meta), path("coverage_plot.png")              , emit: coverage_plot              , optional: true
    tuple val(meta), path("final_mitogenome.annotation.png"), emit: final_mitogenome_annotation, optional: true
    // tuple val(meta), path("final_mitogenome_choice/")       , emit: final_mitogenome_choice    , optional: true
    tuple val(meta), path("final_mitogenome.coverage.png")  , emit: final_mitogenome_coverage  , optional: true
    // tuple val(meta), path("potential_contigs/")             , emit: potential_contigs          , optional: true
    tuple val(meta), path("reads_mapping_and_assembly/")    , emit: reads_mapping_and_assembly , optional: true
    tuple val(meta), path("shared_genes.tsv")               , emit: shared_genes               , optional: true
    tuple val(meta), path("MitoReference")                  , emit: reference_files            , optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.hifiasm.log"), emit: logs 
    tuple val(meta), path("${meta.mt_assembly_prefix}.log"), emit: command_logs 
    path  "versions.yml"                                    , emit: versions
  
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    if (! ["c", "r"].contains(input_mode)) {
        error "r for reads or c for contigs must be specified"
    }

    """
    set +e  # Don't exit on error

    # Run MitoHiFi
    mitohifi.py \\
        -${input_mode} ${hifi_cat} \\
        -f ${ref_fa} \\
        -g ${ref_gb} \\
        -t $task.cpus ${args} \\
        -o ${mito_code}

    exit_code=\$?

    # Rename files
    mv final_mitogenome.fasta ${prefix}.fasta
    mv final_mitogenome.gb ${prefix}.gb
    mv reads_mapping_and_assembly/hifiasm.log ${prefix}.hifiasm.log
    mv .command.log ${prefix}.log
    mv contigs_stats.tsv ${prefix}.contigs_stats.tsv

    mkdir -p MitoReference
    cp $ref_fa $ref_gb MitoReference/
    # Check if the main output was created and rename the fasta header to include the mt_assembly_prefix
    if [ -f "${prefix}.fasta" ]; then
        sed -i "/^>/s/.*/>${prefix}/g" ${prefix}.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    # Create all expected output files for testing
    touch ${prefix}.fasta
    touch ${prefix}.gb
    touch ${prefix}.contigs_stats.tsv
    touch all_potential_contigs.fa
    touch contigs_annotations.png
    touch coverage_plot.png
    touch final_mitogenome.annotation.png
    touch final_mitogenome.coverage.png
    touch shared_genes.tsv

    mkdir contigs_circularization
    mkdir contigs_filtering
    mkdir coverage_mapping
    mkdir final_mitogenome_choice
    mkdir potential_contigs
    mkdir reads_mapping_and_assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//')
    END_VERSIONS
    """
}