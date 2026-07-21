process ASSEMBLY_NO_RESULT {
    tag "$meta.id $reason"
    label 'process_single'

    input:
    tuple val(meta), val(reason)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.assembly_status.tsv"), emit: status
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix ?: meta.id}"
    """
    : > ${prefix}.fasta
    printf 'sample\tstatus\treason\n${prefix}\tdata_limited\t${reason}\n' > ${prefix}.assembly_status.tsv
    printf '"%s":\n    pipeline: "%s"\n' '${task.process}' '${workflow.manifest.version ?: "dev"}' > versions.yml
    """
}
