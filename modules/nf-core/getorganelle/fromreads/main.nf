process GETORGANELLE_FROMREADS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::getorganelle=1.7.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    // `genes` is the sample's curated group gene database (assets/refdb/<group>/
    // <group>_mito_refdb.label.fasta), or [] for a vertebrate / a class with no
    // curated database, in which case no --genes flag is emitted and the command
    // line is byte-identical to the stock-labelling behaviour.
    tuple val(meta), path(fastp), val(organelle_type), path(db), path(genes)

    output:
    tuple val(meta), path("mtdna/${meta.mt_assembly_prefix}.fasta")            , emit: fasta
    tuple val(meta), path("mtdna/${meta.mt_assembly_prefix}.get_org.log.txt")  , emit: log
    path("mtdna/*.selected_graph.gfa")                                               , emit: org_assm_graph,  optional: true
    path("mtdna/*extended_K*.assembly_graph.fastg")                                  , emit: raw_assm_graph,  optional: true
    // extend-* rather than extend-animal_mt: GetOrganelle names these after the
    // LabelDatabase it actually used, which with --genes is the group database
    // (e.g. ...extend-mollusca_mito_refdb.label.fastg). Both are optional outputs,
    // so a stale extend-animal_mt glob would not error -- it would just silently
    // emit nothing for every invertebrate and starve ch_summary_files.
    path("mtdna/*extended_K*.assembly_graph.fastg.extend-*.fastg")           , emit: simp_assm_graph,  optional: true
    path("mtdna/*extended_K*.assembly_graph.fastg.extend-*.csv")             , emit: contig_label,     optional: true
    // path("mtdna/*")                                                         , emit: etc // have only included the files we want above, uncomment if you want everything
    tuple val(meta), path("02_getorganelle_fromreads.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    // --genes overrides the LabelDatabase used to slim the assembly graph; -F still
    // supplies the SeedDatabase that recruits the reads, so read recruitment is
    // unchanged by construction. Mirrors modules/local/getorganelle/reseed.
    def genes_arg = genes ? "--genes ${genes}" : ''
    def effective_args = [args, genes_arg, "--prefix ${prefix}.", "-F ${organelle_type}", "--config-dir ${db}", "-t ${task.cpus}", "-1 ${fastp[0]}", "-2 ${fastp[1]}"].findAll { it?.toString()?.trim() }.join(' ')

    // Use persistent output directory for checkpoint/resume capability
    def checkpoint_base = "${params.outdir}/getorganelle_checkpoints"
    def output_dir = "${checkpoint_base}/${prefix}_getorganelle/mtdna"
    """
    mkdir -p $output_dir

    get_organelle_from_reads.py \\
        $args \\
        $genes_arg \\
        --prefix ${prefix}. \\
        -F $organelle_type \\
        --config-dir $db \\
        -t $task.cpus \\
        -1 ${fastp[0]} \\
        -2 ${fastp[1]} \\
        -o $output_dir

    cat <<-END_TOOL_PARAMS > 02_getorganelle_fromreads.tool_params_mqcrow.html
    <tr><td>GetOrganelle From Reads</td><td><samp>${effective_args}</samp></td><td>Assembles ${organelle_type} from cleaned read pairs for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    wait

    mkdir -p mtdna
    shopt -s nullglob

    # GetOrganelle log is always produced; copy if present, otherwise stub it.
    if [ -f "$output_dir/${meta.mt_assembly_prefix}.get_org.log.txt" ]; then
        cp $output_dir/${meta.mt_assembly_prefix}.get_org.log.txt mtdna/
    else
        touch mtdna/${meta.mt_assembly_prefix}.get_org.log.txt
    fi

    # Assembly graph outputs may be absent if GetOrganelle bailed out very early.
    for f in $output_dir/*.selected_graph.gfa; do cp "\$f" mtdna/; done
    for f in $output_dir/*extended_K*.assembly_graph.fastg; do cp "\$f" mtdna/; done
    for f in $output_dir/*extended_K*.assembly_graph.fastg.extend-*.fastg; do cp "\$f" mtdna/; done
    for f in $output_dir/*extended_K*.assembly_graph.fastg.extend-*.csv; do cp "\$f" mtdna/; done

    # The assembled organelle FASTA (`*1.1.*.fasta`) is only produced when
    # GetOrganelle successfully assembled a contig. When the run finishes
    # without producing one, emit an empty placeholder so downstream stages can
    # detect the "failed to assemble" state instead of crashing on the missing
    # output.
    fasta_files=( $output_dir/${prefix}.*1.1.*.fasta )
    if [ \${#fasta_files[@]} -gt 0 ]; then
        cp "\${fasta_files[@]}" mtdna/
        mv mtdna/${prefix}.*1.1.*.fasta mtdna/${prefix}.fasta
        sed -i "/^>/s/.*/>${prefix}/g" mtdna/${prefix}.fasta
    else
        touch mtdna/${prefix}.fasta
    fi

    wait

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_from_reads.py --version | sed 's/^GetOrganelle v//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    def args   = task.ext.args ?: ''
    def genes_arg = genes ? "--genes ${genes}" : ''
    // The touched extend-animal_mt.* names below still match the extend-* output
    // globs, and deliberately document the no---genes (stock labelling) case.
    def effective_args = [args, genes_arg, "--prefix ${prefix}.", "-F ${organelle_type}", "--config-dir ${db}", "-t ${task.cpus}", "-1 ${fastp[0]}", "-2 ${fastp[1]}"].findAll { it?.toString()?.trim() }.join(' ')
    """
    mkdir -p mtdna
       
    touch "mtdna/${prefix}.fasta"
    touch "mtdna/${prefix}.get_org.log.txt"
    touch "mtdna/${prefix}.selected_graph.gfa"
    touch "mtdna/${prefix}.extended_K21.assembly_graph.fastg"
    touch "mtdna/${prefix}.extended_K21.assembly_graph.fastg.extend-animal_mt.fastg"
    touch "mtdna/${prefix}.extended_K21.assembly_graph.fastg.extend-animal_mt.csv"

    cat <<-END_TOOL_PARAMS > 02_getorganelle_fromreads.tool_params_mqcrow.html
    <tr><td>GetOrganelle From Reads</td><td><samp>${effective_args}</samp></td><td>Assembles ${organelle_type} from cleaned read pairs for ${meta.id}.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """
}
