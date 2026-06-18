// Decide whether a MITOS2 anthozoan annotation is deficient (missing 16S/RNR2 or
// a truncated nad5) and therefore needs the coral fixer. Emits a one-line
// decision file (FIX/PASS) that the annotation subworkflow branches on, so only
// the genuinely broken corals are re-annotated while the correctly annotated
// batch-mates pass through MITOS2 untouched.
process ANNOTATION_QC_GATE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    // The GFF + proteins/ dir MITOS2 produced (MITOS2.out.gff_proteins).
    tuple val(meta), path(gff), path(proteins)

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.coral_qc.txt"), emit: decision
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    annotation_qc_gate.py \\
        --gff ${gff} \\
        --proteins ${proteins} \\
        --out ${meta.mt_assembly_prefix}.coral_qc.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'FIX\\tstub\\n' > ${meta.mt_assembly_prefix}.coral_qc.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
