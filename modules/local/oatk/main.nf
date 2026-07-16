process OATK {
    tag "$meta.id"
    label 'process_medium'

    // Oatk (+ its syncasm assembler and nhmmscan) biocontainer. Pinned via a
    // param so it can be bumped in one place; defaults to the biocontainer build.
    container "${params.oatk_container ?: 'quay.io/biocontainers/oatk:1.0--h577a1d6_1'}"

    input:
    tuple val(meta), path(reads)
    // OatkDB mitochondrial profile-HMM: the <clade>_mito.fam AND its pressed nhmmer
    // index files (.h3f/.h3i/.h3m/.h3p). All are staged together (pass the .fam plus
    // its siblings as a list) because nhmmscan needs the indexes next to the .fam;
    // the .fam itself is resolved from the staged set at runtime.
    path mito_db

    output:
    // Structure-solved MT contigs. optional: Oatk writes nothing when it cannot
    // resolve an organelle contig, so the fallback records a failure rather than
    // crashing the run.
    // The contig is republished as <prefix>.fasta (the same clean name MitoHiFi /
    // GetOrganelle use) so the annotation stage derives the prefix consistently, and
    // the graph as <prefix>.gfa (matching the other assemblers' output names). The
    // native .mito.ctg.fasta is kept for provenance; the .gfa carries the summary's
    // circularity read.
    tuple val(meta), path("${prefix}.fasta"),          emit: fasta, optional: true
    tuple val(meta), path("${prefix}.gfa"),            emit: gfa,   optional: true
    tuple val(meta), path("${prefix}.mito.ctg.fasta"), emit: ctg,   optional: true
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
    # The profile-HMM Oatk reads is the .fam; its .h3* index files must be staged
    # in the same dir (they are, via the module's multi-file db input).
    MITO_DB=\$(ls -1 *.fam 2>/dev/null | head -n1)
    if [ -z "\$MITO_DB" ]; then
        echo "ERROR: no .fam profile-HMM found in the staged Oatk database" >&2
        exit 1
    fi

    oatk \\
        -k ${kmer} \\
        -c ${cov} \\
        -t ${task.cpus} \\
        --nhmmscan \$NHMMSCAN \\
        -m \$MITO_DB \\
        -o ${prefix} \\
        ${args} \\
        ${reads} > ${prefix}.oatk.log 2>&1 || true

    # Republish the structure-solved contig under the clean <prefix>.fasta name that
    # the annotation stage expects. Skipped silently when Oatk resolved no contig
    # (the optional outputs then simply do not exist and the fallback is recorded
    # as a failed assembly downstream).
    if [ -s ${prefix}.mito.ctg.fasta ]; then
        cp ${prefix}.mito.ctg.fasta ${prefix}.fasta
    fi

    # Republish Oatk's native <prefix>.mito.gfa under the clean <prefix>.gfa name the
    # other assemblers use. Skipped silently when no graph was produced.
    if [ -s ${prefix}.mito.gfa ]; then
        cp ${prefix}.mito.gfa ${prefix}.gfa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oatk: \$(oatk --version 2>&1 | head -n1 | sed 's/^oatk //' || echo "unknown")
        nhmmscan: \$(nhmmscan -h 2>&1 | grep -oiE 'HMMER [0-9.]+' | head -n1 | sed 's/HMMER //' || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix ?: meta.id}"
    """
    printf ">%s\\nACGT\\n" ${prefix} > ${prefix}.mito.ctg.fasta
    cp ${prefix}.mito.ctg.fasta ${prefix}.fasta
    printf "S\\tu0\\tACGT\\nL\\tu0\\t+\\tu0\\t+\\t0M\\n" > ${prefix}.gfa
    touch ${prefix}.oatk.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oatk: "stub"
    END_VERSIONS
    """
}
