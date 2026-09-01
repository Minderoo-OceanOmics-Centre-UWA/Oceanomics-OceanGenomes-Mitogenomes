// Compile one run-level held_samples.tsv: a row for every sample that did NOT
// reach ENA submission-ready, with the reason. Two sources feed in as headerless
// per-sample fragments:
//   * pre-QC holds  (proceed_qc=false: <cause>)  -- from EVALUATE_QC_CONDITIONS
//   * table2asn quarantine (FAIL_TABLE2ASN: <blocking_codes>) -- from PARSE_TABLE2ASN_VALIDATION
// Closes the "run reports success with a sixth of the batch quietly missing" gap:
// ch_not_qc_ready and the table2asn quarantine set were console-only .view()s.
//
// Always emits the file (header only when nothing was held), so its presence is a
// reliable end-of-run signal rather than something that appears only on failure.
process COMPILE_HELD_SAMPLES {
    tag 'held_samples'
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path fragments, stageAs: 'fragments/*'

    output:
    path 'held_samples.tsv', emit: held
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    {
        printf 'sample\\tassembly_prefix\\tstage\\treason\\n'
        cat fragments/*.held.tsv 2>/dev/null | sort -u || true
    } > held_samples.tsv

    n=\$(( \$(wc -l < held_samples.tsv) - 1 ))
    echo "[compile_held_samples] \${n} sample(s) held out of submission-ready" >&2

    printf '"%s":\\n    coreutils: "%s"\\n' "${task.process}" "\$(sort --version | head -1 | sed 's/.* //')" > versions.yml
    """

    stub:
    """
    printf 'sample\\tassembly_prefix\\tstage\\treason\\n' > held_samples.tsv
    printf '"%s":\\n    coreutils: "stub"\\n' "${task.process}" > versions.yml
    """
}
