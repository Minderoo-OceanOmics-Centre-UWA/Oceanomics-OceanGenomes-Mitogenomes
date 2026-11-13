// Stage 2 of reference resolution: pick a sample's records out of its curated
// group database by SEQUENCE similarity to its own assembly.
//
// Stage 1 (InvertTaxonGroups.seedDbGroup) narrows a class to one of the
// databases under assets/refdb/. Until this module ran in the reseed, that was
// the only narrowing an invertebrate got: GETORGANELLE_RESEED was handed the
// whole group -- 221 anthozoan genomes, 850 molluscan -- which recruits reads
// from across the phylum, halves effective coverage and shatters the assembly
// graph. Vertebrates never had this problem because findMitoReference resolves a
// single relative and REFERENCE_RANK re-picks it from the reads.
//
// Two modes off the same ranking:
//   'seed'      -> top-n genomes as GETORGANELLE_RESEED's -s, their genes as --genes
//   'reference' -> the single best record as the annotation / QC reference
//
// The reference record is rebuilt from the group's tracked fasta + manifest +
// features by bin/refdb_record.py, not sliced out of a stored .gb -- which is
// what lets every group have a reference rather than only anthozoa, whose .gb
// was the one small enough to track. Always exits 0. Runs in the MITOS2
// BioContainer.
process SELECT_REFERENCE_DB {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // assembly = the mitogenome FASTA to match; refdb_dir = assets/refdb/<group>/.
    tuple val(meta), path(assembly), path(refdb_dir), val(group)
    val mode

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference.gb")        , emit: reference,   optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.seed.fasta")          , emit: seed_fasta,  optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.genedb.fasta")        , emit: label_fasta, optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference_select.txt"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.mt_assembly_prefix}"
    // Seed mode also emits the best record: the reseed needs the top-n as its seed
    // AND that one as the reference the QC checks grade against, and one ranking
    // answers both.
    def out_args = mode == 'seed'
        ? "--out-seed-fasta ${prefix}.seed.fasta --out-label-fasta ${prefix}.genedb.fasta --out-gb ${prefix}.reference.gb"
        : "--out-gb ${prefix}.reference.gb"
    """
    select_reference_db.py \\
        --assembly ${assembly} \\
        --refdb-dir ${refdb_dir} \\
        --group ${group} \\
        ${out_args} \\
        --out-status ${prefix}.reference_select.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>/dev/null | sed -n 's/^blastn: //p')
        python: \$(python --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.mt_assembly_prefix}"
    def stub_out = mode == 'seed'
        ? "printf 'SELECTED_SEED\\t5 of 9 aligned records, 80 label seqs: NC_000000.1[Testidae] best_cov=1.00 best_pid=99.0\\n' > ${prefix}.reference_select.txt\n    touch ${prefix}.seed.fasta ${prefix}.genedb.fasta ${prefix}.reference.gb"
        : "printf 'SELECTED\\tNC_000000.1 Stub coralus [Testidae] cov=1.00 pid=99.0\\n' > ${prefix}.reference_select.txt\n    touch ${prefix}.reference.gb"
    """
    ${stub_out}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
