//
// Subworkflow with functionality specific to the nf-core/oceangenomesmitogenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // 
    // Create channel from input file provided through input
    //

    // Channel
    //     .fromList(samplesheetToList("${input}", "${projectDir}/assets/schema_input.json"))
    //     .map {
    //         meta, fastq_1, fastq_2 ->
    //             if (!fastq_2) {
    //                 return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
    //             } else {
    //                 return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
    //             }
    //     }
    //     .groupTuple()
    //     .map { samplesheet ->
    //         validateInputSamplesheet(samplesheet)
    //     }
    //     .map {
    //         meta, fastqs ->
    //             return [ meta, fastqs.flatten() ]
    //     }
    //     .set { ch_samplesheet }

    // emit:
    // samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    // Capture the workflow metadata handle and params now: inside the onComplete
    // closure the implicit `workflow` and `params` variables resolve to null (they
    // are resolved against the closure delegate, not the script binding), so
    // referencing them there throws a NullPointerException. Capturing them here
    // makes the closure close over valid references.
    def workflow_metadata = workflow
    def trace_report_suffix = params.trace_report_suffix

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

        // errorStrategy 'ignore' keeps the run alive when a process fails, but
        // the default completion summary only reports a count - it never says
        // WHICH samples failed, so silently-dropped tasks (e.g. a reseed whose
        // output was not collected) are easy to miss. Re-surface them loudly
        // here by reading the per-run execution trace, without aborting the run.
        if (workflow_metadata.stats.ignoredCount > 0 || !workflow_metadata.success) {
            def trace_file = file("${outdir}/pipeline_info/execution_trace_${trace_report_suffix}.txt")
            def errored = []
            if (trace_file.exists()) {
                def lines = trace_file.readLines()
                def header = lines ? lines[0].split('\t') : []
                def name_idx = header.findIndexOf { it == 'name' }
                def status_idx = header.findIndexOf { it == 'status' }
                if (name_idx >= 0 && status_idx >= 0) {
                    lines.drop(1).each { line ->
                        def cols = line.split('\t')
                        if (cols.size() > status_idx && cols[status_idx] in ['FAILED', 'ABORTED']) {
                            errored << cols[name_idx]
                        }
                    }
                }
            }
            log.warn "=================================================================="
            log.warn "${workflow_metadata.stats.ignoredCount} process(es) failed and were ignored (run was not aborted)."
            if (errored) {
                log.warn "Failed/ignored tasks:"
                errored.unique().sort().each { log.warn "  - ${it}" }
            }
            log.warn "Review these before trusting downstream results."
            log.warn "Full per-task detail: ${trace_file}"
            log.warn "=================================================================="
        }

        // Memory hints: attempt-scaled processes (OATK, GETORGANELLE_RESEED,
        // GETORGANELLE_FROMREADS, see conf/base.config) request more memory on
        // each retry via task.attempt. task.attempt is scoped to this session,
        // so a later -resume restarts a sample's counter at 1 and re-fails the
        // default low tier before climbing back to the tier that worked last
        // time -- its hash never matches the cached success. Record which
        // tier actually succeeded per sample here, from this run's own trace,
        // so conf/base.config's hintedMemory() can request it directly next
        // time. Purely additive: a missing/unreadable trace is a no-op.
        def hinted_processes = ['OATK', 'GETORGANELLE_RESEED', 'GETORGANELLE_FROMREADS']
        def hints_trace_file = file("${outdir}/pipeline_info/execution_trace_${trace_report_suffix}.txt")
        if (hints_trace_file.exists()) {
            def trace_lines  = hints_trace_file.readLines()
            def trace_header = trace_lines ? trace_lines[0].split('\t') : []
            def name_idx     = trace_header.findIndexOf { it == 'name' }
            def status_idx   = trace_header.findIndexOf { it == 'status' }
            def attempt_idx  = trace_header.findIndexOf { it == 'attempt' }
            def memory_idx   = trace_header.findIndexOf { it == 'memory' }
            if ([name_idx, status_idx, attempt_idx, memory_idx].every { it >= 0 }) {
                def hints_file = file("${outdir}/pipeline_info/memory_hints.json")
                def hints = hints_file.exists() ? new groovy.json.JsonSlurper().parse(hints_file) : [:]
                def hints_updated = false
                trace_lines.drop(1).each { line ->
                    def cols = line.split('\t')
                    if (cols.size() <= memory_idx) { return }
                    if (cols[status_idx] != 'COMPLETED') { return }
                    def attempt = cols[attempt_idx].isInteger() ? cols[attempt_idx].toInteger() : 1
                    if (attempt <= 1) { return }
                    def name_matcher = (cols[name_idx] =~ /:(\w+) \(([^)]+)\)\s*$/)
                    if (!name_matcher.find()) { return }
                    def process_name = name_matcher.group(1)
                    if (!(process_name in hinted_processes)) { return }
                    def sample_id = name_matcher.group(2)
                    def memory_gb = parseTraceMemoryToGb(cols[memory_idx])
                    if (memory_gb == null) { return }
                    if (!hints.containsKey(process_name)) { hints[process_name] = [:] }
                    if (hints[process_name][sample_id] != memory_gb) {
                        hints[process_name][sample_id] = memory_gb
                        hints_updated = true
                    }
                }
                if (hints_updated) {
                    hints_file.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(hints))
                    log.info "Updated memory hints for -resume: ${hints_file}"
                }
            }
        }

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Parse a Nextflow trace `memory` field value (e.g. "200 GB", "6144 MB") into
// whole gigabytes, rounding up so a recorded memory hint never under-provisions
// the next run. Returns null if the value can't be parsed.
//
def parseTraceMemoryToGb(String value) {
    if (!value) { return null }
    def matcher = value.trim() =~ /(?i)^([\d.]+)\s*(B|KB|MB|GB|TB)$/
    if (!matcher.matches()) { return null }
    def amount = matcher.group(1) as Double
    def unit = matcher.group(2).toUpperCase()
    def gb = amount
    if (unit == 'TB')      { gb = amount * 1024 }
    else if (unit == 'MB') { gb = amount / 1024 }
    else if (unit == 'KB') { gb = amount / 1024 / 1024 }
    else if (unit == 'B')  { gb = amount / 1024 / 1024 / 1024 }
    return Math.ceil(gb) as Integer
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() { }

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
