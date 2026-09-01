// Run tRNAscan-SE 2.0 (vertebrate-mitochondrial model) over an EMMA genome FASTA.
//
// tRNAscan-SE's BioContainer is Perl-only (no Python), so this process does just
// the scan and hands its tabular output to TRNA_RESCUE, which parses, guards and
// splices in the stdlib-only psycopg2 container. Only FIX bundles from
// TRNA_RESCUE_GATE reach here.
//
// Fail-safe: errorStrategy 'ignore' plus a `|| true` and an always-created output
// file mean a scan failure yields an empty table, so every target then SKIPs and
// the bundle flows on untouched.
process TRNA_SCAN {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trnascan-se:2.0.12--pl5321h031d066_1' :
        'quay.io/biocontainers/trnascan-se:2.0.12--pl5321h031d066_1' }"

    input:
    // The annotation bundle (EMMA, ND4L/ATP8-rescued); only the genome *.fa is read.
    tuple val(meta), path(annotation, stageAs: 'emma_in/*')

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.trnascan.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def model = task.ext.model ?: 'vert'
    """
    out=${meta.mt_assembly_prefix}.trnascan.tsv
    : > "\$out"

    fa=\$(find emma_in -maxdepth 1 -name '*.fa' | head -n1)
    if [ -n "\$fa" ]; then
        tRNAscan-SE -M ${model} -q -Q ${args} -o "\$out" "\$fa" || : > "\$out"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "tRNAscan-SE": \$(tRNAscan-SE -h 2>&1 | sed -n 's/.*tRNAscan-SE \\([0-9][0-9.]*\\).*/\\1/p' | head -n1)
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.mt_assembly_prefix}.trnascan.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "tRNAscan-SE": "2.0.12"
    END_VERSIONS
    """
}
