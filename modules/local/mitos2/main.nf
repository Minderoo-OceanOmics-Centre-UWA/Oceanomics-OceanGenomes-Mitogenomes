process MITOS2 {
    tag "$meta.id" // tag is used for the publish dir
    label 'process_medium'

    // MITOS2 has no native Singularity build; use the BioContainer for both
    // engines (Setonix pulls the docker:// image via Singularity, like EMMA).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path refdb

    output:
    // Emits the SAME channels as EMMA so the annotation subworkflow can mix()
    // invertebrate (MITOS2) and vertebrate (EMMA) samples without rewiring.
    // Per-gene CDS outputs are optional: a partial/divergent mitogenome may not
    // yield every gene (matches EMMA's behaviour).
    tuple val(meta), path("emma/cds/*CO1*.fa"),  emit: co1_sequences, optional: true
    tuple val(meta), path("emma/cds/*RNR1*.fa"), emit: s12_sequences, optional: true
    tuple val(meta), path("emma/cds/*RNR2*.fa"), emit: s16_sequences, optional: true
    tuple val(meta), path("emma/*"), emit: results
    tuple val(meta), path("07_mitos.tool_params_mqcrow.html"), emit: tool_params
    path "versions_mitos.yml", emit: versions

    script:
        def prefix   = task.ext.prefix ?: meta.mt_assembly_prefix
        def gcode    = task.ext.code ?: (meta.genetic_code ?: 5)
        def refver   = params.mitos_refseq_ver
        def species  = meta.species ?: ''
        // The pinned BioContainer's `mitos` package self-reports version 0.0.0, so
        // derive the version from the container pin instead (bump both together).
        def mitos_version = '2.1.10'
        def mitos_tag     = mitos_version.replaceAll('\\.', '')
        def base_args = (task.ext.args ?: '').toString().trim()
        // GetOrganelle/MitoHiFi assemble a circular molecule; only pass --linear
        // when the assembly is known to be non-circular. For circular (or unknown)
        // topology, run MITOS circular so a gene straddling the linearisation point
        // is annotated as one wrap-around feature instead of being split/dropped.
        def topology_arg = (meta.circular == false) ? '--linear' : ''
        def effective_args = ["runmitos -i ${fasta} -c ${gcode} -r ${refver} -R ${refdb} ${topology_arg} --noplots --best ${base_args}".replaceAll(/ +/, ' ').trim(),
                              "mitos_to_emma.py --bed result.bed --genome ${fasta} --code ${gcode}"].join('; ')

        """
        mkdir -p mitos_raw emma

        # Annotate with MITOS2 (invertebrate mito genetic code by default)
        runmitos \\
            -i ${fasta} \\
            -c ${gcode} \\
            -o mitos_raw \\
            -r ${refver} \\
            -R ${refdb} \\
            ${topology_arg} \\
            --noplots \\
            --best \\
            ${base_args}

        # MITOS writes result.* into the output dir (a per-sequence subdir for
        # multi-contig input); locate the BED regardless of nesting.
        bed=\$(find mitos_raw -name 'result.bed' | head -n1)
        if [ -z "\$bed" ]; then
            echo "ERROR: MITOS2 produced no result.bed under mitos_raw/" >&2
            exit 1
        fi

        mitos_prefix="${prefix}.mitos${mitos_tag}"

        # Adapt MITOS2 output to the EMMA contract (gff + cds/ + proteins/)
        mitos_to_emma.py \\
            --bed "\$bed" \\
            --genome ${fasta} \\
            --prefix "\${mitos_prefix}" \\
            --outdir emma \\
            --code ${gcode} \\
            --species "${species}"

        # Preserve raw MITOS outputs for provenance, but in a subdir so the
        # top-level *.gff glob (used by annotation_stats) only sees the EMMA gff.
        mkdir -p emma/mitos_raw
        cp \$(dirname "\$bed")/result.* emma/mitos_raw/ 2>/dev/null || true

        cat <<-END_TOOL_PARAMS > 07_mitos.tool_params_mqcrow.html
        <tr><td>MITOS2</td><td><samp>${effective_args}</samp></td><td>Annotates the invertebrate mitogenome assembly for ${meta.id} with MITOS2 (genetic code ${gcode}) and reshapes the output to the EMMA contract.</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions_mitos.yml
        "${task.process}":
            mitos: ${mitos_version}
            python: \$(python --version | sed 's/^Python //')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
        def gcode  = task.ext.code ?: (meta.genetic_code ?: 5)
        def base_args = (task.ext.args ?: '').toString().trim()
        def topology_arg = (meta.circular == false) ? '--linear' : ''
        def effective_args = "runmitos -i ${fasta} -c ${gcode} -r ${params.mitos_refseq_ver} -R <refdb> ${topology_arg} --noplots --best ${base_args}".replaceAll(/ +/, ' ').trim()
        """
        mitos_prefix="${prefix}.mitos2110"

        mkdir -p emma/cds emma/proteins emma/mitos_raw

        # Mock EMMA-contract outputs
        touch emma/\${mitos_prefix}.fa
        touch emma/\${mitos_prefix}.gff
        touch emma/cds/CO1.\${mitos_prefix}.fa
        touch emma/cds/RNR1.\${mitos_prefix}.fa
        touch emma/cds/RNR2.\${mitos_prefix}.fa
        touch emma/proteins/MT-CO1.\${mitos_prefix}.fa

        cat <<-END_TOOL_PARAMS > 07_mitos.tool_params_mqcrow.html
        <tr><td>MITOS2</td><td><samp>${effective_args}</samp></td><td>Annotates the invertebrate mitogenome assembly for ${meta.id} with MITOS2 (genetic code ${gcode}) and reshapes the output to the EMMA contract.</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions_mitos.yml
        "${task.process}":
            mitos: "2.1.10"
            python: "3.10.0"
        END_VERSIONS
        """
}
