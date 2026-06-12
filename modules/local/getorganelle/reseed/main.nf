// GetOrganelle "reseed" pass.
//
// A near-copy of modules/nf-core/getorganelle/fromreads that runs a SECOND
// assembly attempt for samples whose first pass failed or produced a poor
// (non-circular / label-database-miss) result because the generic animal_mt
// seed was too divergent from the target. It takes a closely-related
// mitogenome FASTA (downloaded by MITOHIFI_FINDMITOREFERENCE from the resolved
// reference_species_id) and passes it to GetOrganelle as a custom seed (-s),
// while keeping the animal_mt label database (-F animal_mt) for disentangling.
//
// Outputs use a "reseed" suffix (appended to mt_assembly_prefix, no separating
// dot) so they never collide with the first pass in the persistent checkpoint
// directory while keeping the dot-delimited naming convention intact.
process GETORGANELLE_RESEED {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::getorganelle=1.7.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fastp), val(organelle_type), path(db), path(seed)

    output:
    tuple val(meta), path("mtdna/${meta.mt_assembly_prefix}reseed.fasta")            , emit: fasta
    tuple val(meta), path("mtdna/${meta.mt_assembly_prefix}reseed.get_org.log.txt")  , emit: log
    path("mtdna/*.selected_graph.gfa")                                               , emit: org_assm_graph,  optional: true
    path("mtdna/*extended_K*.assembly_graph.fastg")                                  , emit: raw_assm_graph,  optional: true
    path("mtdna/*extended_K*.assembly_graph.fastg.extend-animal_mt.fastg")           , emit: simp_assm_graph, optional: true
    path("mtdna/*extended_K*.assembly_graph.fastg.extend-animal_mt.csv")             , emit: contig_label,    optional: true
    tuple val(meta), path("02b_getorganelle_reseed.tool_params_mqcrow.html")         , emit: tool_params
    path "versions.yml"                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.mt_assembly_prefix}reseed"
    def seed_arg = seed ? "-s ${seed}" : ''
    def effective_args = [args, seed_arg, "--prefix ${prefix}.", "-F ${organelle_type}", "--config-dir ${db}", "-t ${task.cpus}", "-1 ${fastp[0]}", "-2 ${fastp[1]}"].findAll { it?.toString()?.trim() }.join(' ')

    // Use persistent output directory for checkpoint/resume capability
    def checkpoint_base = "${params.outdir}/getorganelle_checkpoints"
    def output_dir = "${checkpoint_base}/${prefix}_getorganelle/mtdna"
    """
    mkdir -p $output_dir

    get_organelle_from_reads.py \\
        $args \\
        $seed_arg \\
        --prefix ${prefix}. \\
        -F $organelle_type \\
        --config-dir $db \\
        -t $task.cpus \\
        -1 ${fastp[0]} \\
        -2 ${fastp[1]} \\
        -o $output_dir

    cat <<-END_TOOL_PARAMS > 02b_getorganelle_reseed.tool_params_mqcrow.html
    <tr><td>GetOrganelle Reseed</td><td><samp>${effective_args}</samp></td><td>Re-assembles ${organelle_type} for ${meta.id} using a closely related seed (${seed}).</td></tr>
    END_TOOL_PARAMS

    wait

    mkdir -p mtdna
    shopt -s nullglob

    # GetOrganelle log is always produced; copy if present, otherwise stub it.
    if [ -f "$output_dir/${prefix}.get_org.log.txt" ]; then
        cp $output_dir/${prefix}.get_org.log.txt mtdna/
    else
        touch mtdna/${prefix}.get_org.log.txt
    fi

    # Assembly graph outputs may be absent if GetOrganelle bailed out very early.
    for f in $output_dir/*.selected_graph.gfa; do cp "\$f" mtdna/; done
    for f in $output_dir/*extended_K*.assembly_graph.fastg; do cp "\$f" mtdna/; done
    for f in $output_dir/*extended_K*.assembly_graph.fastg.extend-animal_mt.fastg; do cp "\$f" mtdna/; done
    for f in $output_dir/*extended_K*.assembly_graph.fastg.extend-animal_mt.csv; do cp "\$f" mtdna/; done

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
    def prefix = "${meta.mt_assembly_prefix}reseed"
    def args   = task.ext.args ?: ''
    def seed_arg = seed ? "-s ${seed}" : ''
    def effective_args = [args, seed_arg, "--prefix ${prefix}.", "-F ${organelle_type}", "--config-dir ${db}", "-t ${task.cpus}", "-1 ${fastp[0]}", "-2 ${fastp[1]}"].findAll { it?.toString()?.trim() }.join(' ')
    """
    mkdir -p mtdna

    touch "mtdna/${prefix}.fasta"
    touch "mtdna/${prefix}.get_org.log.txt"
    touch "mtdna/${prefix}.selected_graph.gfa"
    touch "mtdna/${prefix}.extended_K21.assembly_graph.fastg"
    touch "mtdna/${prefix}.extended_K21.assembly_graph.fastg.extend-animal_mt.fastg"
    touch "mtdna/${prefix}.extended_K21.assembly_graph.fastg.extend-animal_mt.csv"

    cat <<-END_TOOL_PARAMS > 02b_getorganelle_reseed.tool_params_mqcrow.html
    <tr><td>GetOrganelle Reseed</td><td><samp>${effective_args}</samp></td><td>Re-assembles ${organelle_type} for ${meta.id} using a closely related seed (${seed}).</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """
}
