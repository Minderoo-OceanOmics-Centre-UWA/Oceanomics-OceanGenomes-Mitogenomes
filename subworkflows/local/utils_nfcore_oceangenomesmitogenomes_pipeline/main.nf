//
// Subworkflow with functionality specific to the OceanOmics OceanGenomes Mitogenomes pipeline
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
// One entry per tool: [enabled, in-text citation, bibliography entry]. Kept as a
// single list so the citation sentence and the reference list can't drift apart.
// `enabled` mirrors the params that gate each branch, so a run only cites the
// tools it could actually have used. EMMA (vertebrates) and MITOS2 (inverts) are
// routed per sample on meta.invertebrates, so both are cited whenever
// annotation is on.
def toolReferences() {
    def assembly_getorg = !params.skip_mitogenome_assembly_getorg
    def assembly_hifi   = !params.skip_mitogenome_assembly_hifi
    def annotation      = !params.skip_mitogenome_annotation
    def ena_prep        = !params.skip_ena_submission_prep
    return [
        [assembly_getorg, "fastp (Chen et al. 2018)",
            "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. doi: <a href=\"https://doi.org/10.1093/bioinformatics/bty560\">10.1093/bioinformatics/bty560</a></li>"],
        [assembly_getorg, "GetOrganelle (Jin et al. 2020)",
            "<li>Jin, J.-J., Yu, W.-B., Yang, J.-B., Song, Y., dePamphilis, C. W., Yi, T.-S., & Li, D.-Z. (2020). GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology, 21, 241. doi: <a href=\"https://doi.org/10.1186/s13059-020-02154-5\">10.1186/s13059-020-02154-5</a></li>"],
        [assembly_hifi, "MitoHiFi (Uliano-Silva et al. 2023)",
            "<li>Uliano-Silva, M., Ferreira, J. G. R. N., Krasheninnikova, K., Darwin Tree of Life Consortium, Formenti, G., Abueg, L., Torrance, J., Myers, E. W., Durbin, R., Blaxter, M., & McCarthy, S. A. (2023). MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio high fidelity reads. BMC Bioinformatics, 24, 288. doi: <a href=\"https://doi.org/10.1186/s12859-023-05385-y\">10.1186/s12859-023-05385-y</a></li>"],
        [assembly_hifi && params.enable_oatk_fallback, "oatk (Zhou et al. 2025)",
            "<li>Zhou, C., Brown, M., Blaxter, M., Darwin Tree of Life Consortium, McCarthy, S. A., & Durbin, R. (2025). Oatk: a de novo assembly tool for complex plant organelle genomes. Genome Biology, 26, 235. doi: <a href=\"https://doi.org/10.1186/s13059-025-03676-6\">10.1186/s13059-025-03676-6</a></li>"],
        [!params.skip_mitogenome_depth || assembly_hifi, "minimap2 (Li 2018)",
            "<li>Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100. doi: <a href=\"https://doi.org/10.1093/bioinformatics/bty191\">10.1093/bioinformatics/bty191</a></li>"],
        [assembly_hifi, "SAMtools (Danecek et al. 2021)",
            "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: <a href=\"https://doi.org/10.1093/gigascience/giab008\">10.1093/gigascience/giab008</a></li>"],
        [annotation, "EMMA (Al Arab et al. 2017)",
            "<li>Al Arab, M., Höner zu Siederdissen, C., Tout, K., Sahyoun, A. H., Stadler, P. F., & Bernt, M. (2017). Accurate annotation of protein-coding genes in mitochondrial genomes. Molecular Phylogenetics and Evolution, 106, 209-216. doi: <a href=\"https://doi.org/10.1016/j.ympev.2016.09.024\">10.1016/j.ympev.2016.09.024</a></li>"],
        [annotation, "MITOS2 (Donath et al. 2019)",
            "<li>Donath, A., Jühling, F., Al-Arab, M., Bernhart, S. H., Reinhardt, F., Stadler, P. F., Middendorf, M., & Bernt, M. (2019). Improved annotation of protein-coding genes boundaries in metazoan mitochondrial genomes. Nucleic Acids Research, 47(20), 10543-10552. doi: <a href=\"https://doi.org/10.1093/nar/gkz833\">10.1093/nar/gkz833</a></li>"],
        [annotation, "tRNAscan-SE (Chan et al. 2021)",
            "<li>Chan, P. P., Lin, B. Y., Mak, A. J., & Lowe, T. M. (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Research, 49(16), 9077-9096. doi: <a href=\"https://doi.org/10.1093/nar/gkab688\">10.1093/nar/gkab688</a></li>"],
        [annotation, "BLAST+ (Camacho et al. 2009)",
            "<li>Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. doi: <a href=\"https://doi.org/10.1186/1471-2105-10-421\">10.1186/1471-2105-10-421</a></li>"],
        [annotation, "the NCBI Taxonomy database (Schoch et al. 2020)",
            "<li>Schoch, C. L., et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. Database, 2020, baaa062. doi: <a href=\"https://doi.org/10.1093/database/baaa062\">10.1093/database/baaa062</a></li>"],
        [true, "Biopython (Cock et al. 2009)",
            "<li>Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btp163\">10.1093/bioinformatics/btp163</a></li>"],
        [ena_prep, "EMBOSS seqret (Rice et al. 2000)",
            "<li>Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: the European Molecular Biology Open Software Suite. Trends in Genetics, 16(6), 276-277. doi: <a href=\"https://doi.org/10.1016/S0168-9525(00)02024-2\">10.1016/S0168-9525(00)02024-2</a></li>"],
        [ena_prep, "ENA Webin-CLI",
            "<li>European Nucleotide Archive. Webin-CLI. <a href=\"https://github.com/enasequence/webin-cli\">https://github.com/enasequence/webin-cli</a></li>"],
        [true, "MultiQC (Ewels et al. 2016)",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btw354\">10.1093/bioinformatics/btw354</a></li>"],
    ].findAll { it[0] }
}

def toolCitationText() {
    return "Tools used in the workflow included: " + toolReferences().collect { it[1] }.join(', ') + "."
}

def toolBibliographyText() {
    return toolReferences().collect { it[2] }.join(' ')
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
    meta["tool_citations"] = toolCitationText()
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
