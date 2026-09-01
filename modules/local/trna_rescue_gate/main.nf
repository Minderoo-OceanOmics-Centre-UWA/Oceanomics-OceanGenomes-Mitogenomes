// Decide whether an EMMA annotation is a candidate for the tRNA rescue.
// EMMA's covariance model periodically misses a tRNA that is physically present
// in the assembly on an otherwise complete, correctly ordered vertebrate
// mitogenome. Emits a one-line decision file (FIX\t<targets> | PASS\t-) that the
// annotation subworkflow branches on, so only the recoverable cases go to
// TRNA_RESCUE and everything else passes through untouched.
process TRNA_RESCUE_GATE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    // The annotation bundle (EMMA, ND4L/ATP8-rescued); only the *.gff is read.
    tuple val(meta), path(annotation, stageAs: 'emma_in/*')

    output:
    tuple val(meta), path("${meta.mt_assembly_prefix}.trna_rescue_qc.txt"), emit: decision
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gff=\$(find emma_in -maxdepth 1 -name '*.gff' | head -n1)
    if [ -z "\$gff" ]; then
        printf 'PASS\\t-\\n' > ${meta.mt_assembly_prefix}.trna_rescue_qc.txt
    else
        trna_rescue_gate.py --gff "\$gff" --out ${meta.mt_assembly_prefix}.trna_rescue_qc.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'PASS\\t-\\n' > ${meta.mt_assembly_prefix}.trna_rescue_qc.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10.0"
    END_VERSIONS
    """
}
