// Harness: a database upload that fails must fail its task, without stopping the run, and be
// listed in sql_uploaded_data/FAILED_UPLOADS.tsv by the completion handler.
//
// It must run as a script in the repository root, so bin/, nextflow.config and
// conf/base.config's errorStrategy are the pipeline's own (-main-script only applies to
// pulled projects). From the repository root:
//   sed -e "s#'../../#'./#g" tests/upload_failure/main.nf > upload_failure_harness.nf
//   nextflow run upload_failure_harness.nf -c pawsey_profile.config -profile singularity \
//       --outdir <dir>
// dead.cfg points at a closed port, so nothing is written to any database.
// Run 2026-11-04 on Setonix: the upload failed (exit 1), was retried once and ignored,
// STILL_RUNS completed, the run ended "completed successfully, but with errored process(es)"
// (exit 0), FAILED_UPLOADS.tsv named the task and its work dir, and -resume re-ran only it.
//
// Expected: PUSH_MTDNA_ASSM_RESULTS fails twice (one retry), is ignored; STILL_RUNS completes;
// the run ends successfully; <outdir>/sql_uploaded_data/FAILED_UPLOADS.tsv names the upload.

nextflow.enable.dsl = 2

include { PUSH_MTDNA_ASSM_RESULTS } from '../../modules/local/upload_results/mtdna/main'
include { PIPELINE_COMPLETION     } from '../../subworkflows/local/utils_nfcore_oceangenomesmitogenomes_pipeline'

// Stands in for the rest of the pipeline: it must still run when the upload fails.
process STILL_RUNS {
    label 'process_single'
    output:
    path 'still_ran.txt'
    script:
    """
    echo ok > still_ran.txt
    """
}

workflow {
    def here = "${projectDir}/tests/upload_failure"
    def meta = [id: 'OG1', mt_assembly_prefix: 'OG1.ilmn.250101.getorg1770', circular: true]
    PUSH_MTDNA_ASSM_RESULTS(
        Channel.of([meta,
                    file("${here}/OG1.ilmn.250101.getorg1770.fasta"),
                    file("${here}/OG1.ilmn.250101.getorg1770.get_org.log.txt"),
                    file("${projectDir}/assets/placeholders/empty_mito_depth.tsv")]),
        file("${here}/dead.cfg"))
    STILL_RUNS()
    PIPELINE_COMPLETION(params.email, params.email_on_fail, params.plaintext_email, params.outdir,
                        params.monochrome_logs, params.hook_url, Channel.empty())
}
