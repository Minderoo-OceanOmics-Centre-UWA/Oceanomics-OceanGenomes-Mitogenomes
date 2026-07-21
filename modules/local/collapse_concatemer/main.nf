process COLLAPSE_CONCATEMER {
    tag "$meta.id"
    label 'process_low'

    // Reuse the pinned MitoHiFi image: it ships makeblastdb / blastn / python3,
    // which is everything collapse_concatemer.py needs to re-verify and collapse
    // a head-to-tail multimer. No extra container required.
    container "${params.mitohifi_container}"

    input:
    // evidence is the <prefix>.circularity_check.tsv emitted by the circularity
    // check; it carries anomaly_type + suggested_trim_region.
    tuple val(meta), path(fasta), path(evidence)

    output:
    // Collapsed FASTA is a drop-in replacement for the assembly FASTA (a genuine
    // concatemer is rewritten to its monomer; everything else is passed through
    // byte-for-byte), written under collapsed/ to avoid clashing with the staged
    // input of the same name.
    tuple val(meta), path("collapsed/${prefix}.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.concatemer_collapse.tsv"), emit: report
    tuple val(meta), path("${prefix}.post_curation_check.tsv"), emit: evidence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix ?: meta.id}"
    """
    mkdir -p collapsed

    collapse_concatemer.py \\
        --fasta ${fasta} \\
        --evidence ${evidence} \\
        --sample ${prefix} \\
        --out-fasta collapsed/${prefix}.fasta \\
        --out-report ${prefix}.concatemer_collapse.tsv \\
        --out-evidence ${prefix}.post_curation_check.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse_concatemer: "1.0.0"
        python: \$(python --version 2>&1 | sed 's/Python //')
        blastn: \$(blastn -version 2>&1 | head -n1 | sed 's/blastn: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix ?: meta.id}"
    """
    mkdir -p collapsed
    cp ${fasta} collapsed/${prefix}.fasta
    printf "sample\\taction\\toriginal_length\\tcollapsed_length\\treference_length\\tinferred_monomer_period\\ttail_identity\\treason\\n${prefix}\\tpassthrough\\tNA\\tNA\\tNA\\tNA\\tNA\\tstub\\n" > ${prefix}.concatemer_collapse.tsv
    cp ${evidence} ${prefix}.post_curation_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse_concatemer: "1.0.0"
    END_VERSIONS
    """
}
