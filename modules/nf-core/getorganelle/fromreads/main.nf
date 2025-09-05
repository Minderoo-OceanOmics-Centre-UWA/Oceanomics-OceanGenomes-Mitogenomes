process GETORGANELLE_FROMREADS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::getorganelle=1.7.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fastp), val(organelle_type), path(db)

    output:
    tuple val(meta), path("mtdna/${meta.mt_assembly_prefix}.fasta")            , emit: fasta
    tuple val(meta), path("mtdna/${meta.mt_assembly_prefix}.get_org.log.txt")  , emit: log
    path("mtdna/*")                                                         , emit: etc // could alter this to just move the files we want, but then you need to go back into workdir for others
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"

    """
    get_organelle_from_reads.py \\
        $args \\
        --prefix ${prefix}. \\
        -F $organelle_type \\
        --config-dir $db \\
        -t $task.cpus \\
        -1 ${fastp[0]} \\
        -2 ${fastp[1]} \\
        -o mtdna

    wait
    
    # Move and rename the output file to include the version in the filename
    mv mtdna/${prefix}.*1.1.*.fasta mtdna/${prefix}.fasta
    sed -i "/^>/s/.*/>${prefix}/g" mtdna/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_from_reads.py --version | sed 's/^GetOrganelle v//g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def base_prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    """
    mkdir -p mtdna
       
    touch "mtdna/${prefix}.fasta"
    touch "mtdna/stub_output.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """
}