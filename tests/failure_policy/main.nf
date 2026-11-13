nextflow.enable.dsl = 2

process SAMPLE_TASK {
    tag "$sample"

    input:
    tuple val(sample), val(mode)

    output:
    tuple val(sample), path("${sample}.txt")

    script:
    // 'lost_ack' models a dropped sbatch acknowledgement: the job really ran and
    // succeeded, but Nextflow saw sbatch's own exit 1. It is exit 1 on the first
    // attempt only, which is indistinguishable by exit code from
    // 'deterministic_failure' -- that is exactly why the policy grants one retry
    // unconditionally rather than classifying by an exit-code allowlist.
    def exit_code = mode == 'deterministic_failure' ? 1 :
        mode == 'lost_ack' && task.attempt == 1 ? 1 :
        mode == 'transient_then_success' && task.attempt == 1 ? 137 :
        mode == 'transient_failure' ? 137 : 0
    """
    if [ ${exit_code} -ne 0 ]; then
        exit ${exit_code}
    fi
    echo '${sample}' > ${sample}.txt
    """
}

process PREPARE_LCA_DATABASES {
    input:
    val mode

    output:
    path 'shared.txt'

    script:
    def exit_code = mode == 'failure' ? 1 : 0
    """
    if [ ${exit_code} -ne 0 ]; then
        exit ${exit_code}
    fi
    touch shared.txt
    """
}

workflow SAMPLE_POLICY {
    take:
    samples

    main:
    SAMPLE_TASK(samples)

    emit:
    results = SAMPLE_TASK.out
}

workflow SHARED_POLICY {
    take:
    mode

    main:
    PREPARE_LCA_DATABASES(mode)

    emit:
    result = PREPARE_LCA_DATABASES.out
}

// Direct entry point for validating the policy without nf-test.
workflow {
    if (params.scenario == 'shared_failure') {
        SHARED_POLICY(Channel.value('failure'))
    } else if (params.scenario == 'lost_ack') {
        SAMPLE_POLICY(Channel.of(['ack', 'lost_ack']))
    } else if (params.scenario == 'transient_success') {
        SAMPLE_POLICY(Channel.of(['retry', 'transient_then_success']))
    } else if (params.scenario == 'transient_exhausted') {
        SAMPLE_POLICY(Channel.of(
            ['good', 'success'],
            ['exhausted', 'transient_failure']
        ))
    } else {
        SAMPLE_POLICY(Channel.of(
            ['good', 'success'],
            ['bad', 'deterministic_failure']
        ))
    }
}
