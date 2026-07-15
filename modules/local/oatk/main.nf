process OATK {
    tag "$meta.id"
    label 'process_medium'

    // Oatk (+ its syncasm assembler and nhmmscan) biocontainer. Pinned via a
    // param so it can be bumped in one place; defaults to the biocontainer build.
    container "${params.oatk_container ?: 'quay.io/biocontainers/oatk:1.0--h577a1d6_1'}"

    input:
    tuple val(meta), path(reads)
    // OatkDB mitochondrial profile-HMM (.fam) for the sample's clade, plus its
    // pressed nhmmer index files (.h3f/.h3i/.h3m/.h3p) staged alongside it.
    path mito_db

    output:
    // Structure-solved MT contigs. optional: Oatk writes nothing when it cannot
    // resolve an organelle contig, so the fallback records a failure rather than
    // crashing the run.
    tuple val(meta), path("${prefix}.mito.ctg.fasta"), emit: fasta, optional: true
    tuple val(meta), path("${prefix}.mito.gfa"),       emit: gfa,   optional: true
    tuple val(meta), path("${prefix}.oatk.log"),       emit: log
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix ?: meta.id}"
    // Syncmer size (HiFi default 1001) and coverage threshold (~5-10x nuclear
    // coverage per the Oatk docs). Both overridable via params without editing
    // this module.
    def kmer = params.oatk_syncmer_size ?: 1001
    def cov  = params.oatk_syncmer_coverage ?: 30
    """
    # nhmmscan ships in the Oatk container; resolve its path for the --nhmmscan arg.
    NHMMSCAN=\$(command -v nhmmscan)

    oatk \\
        -k ${kmer} \\
        -c ${cov} \\
        -t ${task.cpus} \\
        --nhmmscan \$NHMMSCAN \\
        -m ${mito_db} \\
        -o ${prefix} \\
        ${args} \\
        ${reads} > ${prefix}.oatk.log 2>&1 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oatk: \$(oatk --version 2>&1 | head -n1 | sed 's/^oatk //' || echo "unknown")
        nhmmscan: \$(nhmmscan -h 2>&1 | grep -oiE 'HMMER [0-9.]+' | head -n1 | sed 's/HMMER //' || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix ?: meta.id}"
    """
    touch ${prefix}.mito.ctg.fasta ${prefix}.mito.gfa ${prefix}.oatk.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oatk: "stub"
    END_VERSIONS
    """
}
