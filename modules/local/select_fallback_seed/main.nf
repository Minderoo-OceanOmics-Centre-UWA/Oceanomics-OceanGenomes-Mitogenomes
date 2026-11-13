// Select a bounded, read-supported seed panel for an invertebrate whose first
// GetOrganelle pass produced no sequence to use for assembly-based ranking.
process SELECT_FALLBACK_SEED {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(refdb_dir), val(group)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference.gb")         , emit: reference,   optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.seed.fasta")           , emit: seed_fasta,  optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.genedb.fasta")         , emit: label_fasta, optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.fallback_select.txt")  , emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.mt_assembly_prefix
    def readArgs = reads instanceof List ? reads.join(' ') : reads
    def taxon = meta.nominal_species_id ?: meta.species_id ?: ''
    def family = meta.family ?: ''
    def order = meta.order ?: ''
    def taxClass = meta.class ?: ''
    """
    select_fallback_seed.py \
        --reads ${readArgs} \
        --refdb-dir ${refdb_dir} \
        --group ${group} \
        --taxon '${taxon}' \
        --family '${family}' \
        --order '${order}' \
        --class-name '${taxClass}' \
        --out-seed-fasta ${prefix}.seed.fasta \
        --out-label-fasta ${prefix}.genedb.fasta \
        --out-gb ${prefix}.reference.gb \
        --out-status ${prefix}.fallback_select.txt \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>/dev/null | sed -n 's/^blastn: //p')
        python: \$(python --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.mt_assembly_prefix
    """
    touch ${prefix}.seed.fasta ${prefix}.genedb.fasta ${prefix}.reference.gb
    printf 'SELECTED_SEED\tfallback=taxonomy_reads tier=order candidates=5 sampled_reads=100 chosen=NC_000001.1\n' > ${prefix}.fallback_select.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
