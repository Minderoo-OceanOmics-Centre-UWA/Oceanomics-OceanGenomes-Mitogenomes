// Flag when the mitogenome reference chosen for a sample is not relevant to its
// assembly. The reference is resolved by MITOHIFI_FINDMITOREFERENCE from the
// sample's species *label*, so a wrong/coarse label yields a wrong-family
// reference that silently degrades seeding + the coral annotation fix. This
// module BLASTs the reference against the assembly and writes a one-line
// PASS/MISMATCH/UNKNOWN flag (label- and taxonomy-DB-free). Always exits 0 so a
// bad reference only records a review flag, never breaks the run. Runs in the
// MITOS2 BioContainer (provides blastn + biopython).
process REFERENCE_RELEVANCE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // assembly = the published mitogenome FASTA, reference_gb = the resolved ref.
    tuple val(meta), path(assembly), path(reference_gb)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference_relevance.txt"), emit: flag
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    reference_relevance_check.py \\
        --assembly ${assembly} \\
        --reference-gb ${reference_gb} \\
        --out ${meta.mt_assembly_prefix}.reference_relevance.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>/dev/null | sed -n 's/^blastn: //p')
        python: \$(python --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'PASS\\tstub\\n' > ${meta.mt_assembly_prefix}.reference_relevance.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
