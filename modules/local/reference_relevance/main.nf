// Grade how well the mitogenome reference chosen for a sample corresponds to its
// assembly. The reference is resolved by MITOHIFI_FINDMITOREFERENCE from the
// sample's species *label*, so a wrong/coarse label yields a wrong-family
// reference that silently degrades seeding + the coral annotation fix. This
// module BLASTs the reference against the assembly and writes a one-line
// PASS/DIVERGENT/MISMATCH/UNKNOWN flag (label- and taxonomy-DB-free):
//   DIVERGENT = right molecule, distant relative (advisory; a closer reference
//               would seed and annotate better)
//   MISMATCH  = the reference neither covers nor matches, i.e. the wrong reference
// Always exits 0 so a bad reference only records a review flag, never breaks the
// run. Runs in the MITOS2 BioContainer (provides blastn + biopython).
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
    // Sample species drives the congeneric veto: a same-genus reference is the best
    // obtainable, so it is never called MISMATCH however low the identity runs.
    def species = (meta.nominal_species_id ?: meta.reference_species_id ?: '').toString().trim()
    def species_arg = species ? "--sample-species '${species}'" : ''
    // Taxon-aware identity floor. The 88.0 default was calibrated on corals, whose
    // mtDNA evolves far more slowly than vertebrate mtDNA; teleost *congeners*
    // routinely align at 78-88%, so keeping 88.0 for fish flags good assemblies.
    // Only Cnidaria keeps that validated coral number -- Porifera has no evidence
    // of sharing Cnidaria's unusually slow substitution rate, and other invert
    // phyla (Mollusca, Arthropoda, Echinodermata) default to the vertebrate floor
    // as an untuned starting point; bilaterian invertebrate mtDNA (especially
    // arthropod/mollusc) often evolves *faster* than vertebrate mtDNA, so 82.0 may
    // still be too strict -- revisit once real divergence data from this batch is
    // available.
    def min_pid = InvertTaxonGroups.isCoralFixEligible(meta.class) ? 88.0 : 82.0
    """
    reference_relevance_check.py \\
        --assembly ${assembly} \\
        --reference-gb ${reference_gb} \\
        ${species_arg} \\
        --min-pid ${min_pid} \\
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
