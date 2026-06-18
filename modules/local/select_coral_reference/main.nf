// Pick the per-sample coral annotation reference from the curated Anthozoa
// mitogenome DB by SEQUENCE similarity to the assembly, instead of from the
// species label (findMitoReference by name). Label-free, so a wrong/coarse label
// can no longer hand a wrong-family reference to the coral annotation fixer.
// Emits the chosen record's GenBank (with the 16S + nad5 features coral_fix needs)
// plus a status line. Always exits 0. Runs in the MITOS2 BioContainer.
process SELECT_CORAL_REFERENCE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitos:2.1.10--pyhdfd78af_0' :
        'quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0' }"

    input:
    // assembly = the mitogenome FASTA to match; db_gb = assets/coral_mito_refdb.gb.
    tuple val(meta), path(assembly)
    path db_gb

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference.gb"),         emit: reference, optional: true
    tuple val(meta), path("${meta.mt_assembly_prefix}.reference_select.txt"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    select_coral_reference.py \\
        --assembly ${assembly} \\
        --db-gb ${db_gb} \\
        --out-gb ${meta.mt_assembly_prefix}.reference.gb \\
        --out-status ${meta.mt_assembly_prefix}.reference_select.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>/dev/null | sed -n 's/^blastn: //p')
        python: \$(python --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'SELECTED\\tNC_000000.1 Stub coralus [Testidae] cov=1.00 pid=99.0\\n' > ${meta.mt_assembly_prefix}.reference_select.txt
    touch ${meta.mt_assembly_prefix}.reference.gb
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
