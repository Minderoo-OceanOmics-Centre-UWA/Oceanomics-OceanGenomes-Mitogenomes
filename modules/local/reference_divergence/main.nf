// Flag, BEFORE assembly, when the mitogenome reference is not a close relative of
// the sample. MITOHIFI_FINDMITOREFERENCE walks the sample's NCBI lineage and
// downloads the first complete mitogenome it finds, so a species with no
// congeneric record (common for deep-sea taxa) silently gets a divergent
// reference. MitoHiFi then recruits reads against it and drops the most divergent
// gene blocks, yielding a clean-looking but gene-incomplete collapse. This module
// compares taxonomy only (genus, and family when available) and writes a one-line
// CONGENERIC|CONFAMILIAL|DIFFERENT_FAMILY|NON_CONGENERIC|UNKNOWN flag. It is the
// pre-assembly, taxonomy-based complement to the post-assembly REFERENCE_RELEVANCE
// BLAST check. Always exits 0 so a bad reference only records a review flag, never
// breaks the run. Runs in the MITOS2 BioContainer (provides biopython).
process REFERENCE_DIVERGENCE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // reference_gb = the resolved per-sample reference GenBank (taxonomy intact).
    tuple val(meta), path(reference_gb)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference_divergence.txt"), emit: flag
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // nominal species drives the genus comparison; class enriches the message;
    // family (absent from current samplesheets) upgrades the tier when present.
    def species = (meta.nominal_species_id ?: meta.reference_species_id ?: '').toString().trim()
    def klass = (meta.class ?: '').toString().trim()
    def family = (meta.family ?: '').toString().trim()
    def taxon_order = (meta.order ?: '').toString().trim()
    def class_arg = klass ? "--sample-class '${klass}'" : ''
    def family_arg = family ? "--sample-family '${family}'" : ''
    // Enables the CROSS_ORDER tier when the samplesheet supplies an order; absent
    // on older samplesheets, in which case grading falls back to genus/family.
    def order_arg = taxon_order ? "--sample-order '${taxon_order}'" : ''
    """
    reference_divergence_check.py \\
        --reference-gb ${reference_gb} \\
        --sample-species '${species}' \\
        ${class_arg} \\
        ${family_arg} \\
        ${order_arg} \\
        --out ${meta.mt_assembly_prefix}.reference_divergence.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'CONGENERIC\\tstub\\n' > ${meta.mt_assembly_prefix}.reference_divergence.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
