// Compile one run-level held_samples.tsv: a row for every sample that did NOT
// reach ENA submission-ready, with the reason. Three sources feed in as headerless
// per-sample fragments:
//   * pre-annotation drops (assembly_failed:empty_fasta, under_min_length:<bp><<floor>)
//     -- from the annotation-input filters in the main workflow
//   * pre-QC holds  (proceed_qc=false: <cause>)  -- from EVALUATE_QC_CONDITIONS
//   * table2asn quarantine (FAIL_TABLE2ASN: <blocking_codes>) -- from PARSE_TABLE2ASN_VALIDATION
// Closes the "run reports success with a sixth of the batch quietly missing" gap:
// ch_not_qc_ready and the table2asn quarantine set were console-only .view()s.
//
// The grain is one row per ASSEMBLY PER STAGE, not one row per sample, and the
// fragment filenames carry the stage for exactly that reason. A sample rescued by
// reseed shows a PRE_ANNOTATION row for its failed first pass alongside its later
// row, which is correct at this file's grain but means "count the rows" is not
// "count held samples". The `stage` column is what separates them.
//
// The completeness check is the part that does not depend on knowing where a sample
// left. Three fragment sources fix three known drops; enumerating every point a
// sample can leave is not tractable, because any ignored task failure anywhere in the
// assembly subworkflow produces the same shape -- a sample with no assembly in any
// filtered channel, no QC row, and no DB row. Comparing the samplesheet against
// held + submission-ready needs to know nothing about the cause.
//
// It WARNS rather than aborts, following this codebase's convention for this
// situation (see the diagnostic view in attachCircularityEvidence): a reporting gap
// discovered at the very end of a multi-day run should not throw the run away, it
// should be impossible to miss in the log.
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
    path samplesheet_ogs
    path ena_run_summary

    output:
    path 'held_samples.tsv', emit: held
    path 'run_completeness.txt', emit: completeness
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
    echo "[compile_held_samples] \${n} row(s) held out of submission-ready" >&2

    # Completeness: samplesheet - (held + submission-ready). Column 1 of
    # held_samples.tsv is the sample id; og_id is column 2 of ena_run_summary.tsv and
    # submission_ready is column 11, both named in the header rather than assumed by
    # position, so a column added upstream cannot silently shift the check.
    compile_held_completeness.py \\
        --samplesheet-ogs ${samplesheet_ogs} \\
        --held held_samples.tsv \\
        --ena-run-summary ${ena_run_summary} \\
        --output run_completeness.txt

    printf '"%s":\\n    coreutils: "%s"\\n    python: "%s"\\n' \\
        "${task.process}" \\
        "\$(sort --version | head -1 | sed 's/.* //')" \\
        "\$(python --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    printf 'sample\\tassembly_prefix\\tstage\\treason\\n' > held_samples.tsv
    printf 'all %s samplesheet sample(s) accounted for\\n' 0 > run_completeness.txt
    printf '"%s":\\n    coreutils: "stub"\\n    python: "stub"\\n' "${task.process}" > versions.yml
    """
}
