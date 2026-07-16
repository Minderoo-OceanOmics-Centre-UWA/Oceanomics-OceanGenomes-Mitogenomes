#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENA_FLATFILE       } from './modules/local/genome_qc/ena_flatfile'
include { ENA_EMBL_PREFLIGHT } from './modules/local/genome_qc/ena_embl_preflight'
include { WEBIN_VALIDATE     } from './modules/local/genome_qc/webin_validate'
include { ENA_VALIDATION_RESULT  } from './modules/local/genome_qc/ena_validation_result'
include { ENA_VALIDATION_SUMMARY } from './modules/local/genome_qc/ena_validation_summary'
include { UPLOAD_ENA_RESULTS     } from './subworkflows/local/upload_results_mito'
include { MULTIQC                } from './modules/nf-core/multiqc/main'

def requiredValue(row, String column, rowLabel) {
    def value = row[column]?.toString()?.trim()
    if (!value) {
        error "ENA input ${rowLabel} is missing required column '${column}'."
    }
    return value
}

def safeIdentifier(String value, String column, rowLabel) {
    if (!(value ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "ENA input ${rowLabel} has invalid ${column} '${value}'. Use letters, numbers, dots, underscores, or hyphens."
    }
    return value
}

def resolveInputPath(inputSheet, String rawPath, String column, rowLabel) {
    def candidate = java.nio.file.Paths.get(rawPath)
    def resolved = candidate.isAbsolute() ? candidate.normalize() : inputSheet.parent.resolve(candidate).normalize()
    if (!java.nio.file.Files.isRegularFile(resolved)) {
        error "ENA input ${rowLabel} ${column} file does not exist: ${resolved}"
    }
    return file(resolved.toString(), checkIfExists: true)
}

def readNamedStatus(statusFile) {
    def rows = statusFile.readLines()
        .findAll { line -> line?.trim() }
        .collect { line -> line.split('\\t', -1) as List }
    if (rows.size() < 2) {
        return [status: 'MALFORMED_STATUS', reason: 'missing_header_or_data']
    }
    def header = rows[0]
    def statusIndex = header.indexOf('status')
    def reasonIndex = header.indexOf('reason')
    if (statusIndex < 0 || statusIndex >= rows[1].size()) {
        return [status: 'MALFORMED_STATUS', reason: 'missing_status_column']
    }
    def status = rows[1][statusIndex]?.trim()
    def reason = reasonIndex >= 0 && reasonIndex < rows[1].size() ? rows[1][reasonIndex]?.trim() : ''
    return [status: status ?: 'MALFORMED_STATUS', reason: reason ?: 'unspecified']
}

def table2asnGate(statusFile) {
    def parsed = readNamedStatus(statusFile)
    if (parsed.status == 'PASS') {
        return [status: 'PASS', reason: 'table2asn_pass']
    }
    if (parsed.status == 'MALFORMED_STATUS') {
        return [status: 'SKIP_TABLE2ASN', reason: parsed.reason]
    }
    return [status: 'SKIP_TABLE2ASN', reason: "table2asn_${parsed.status}".toLowerCase()]
}

workflow {
    if (!(params.ena_mode in ['convert_validate', 'validate'])) {
        error "--ena_mode must be either 'convert_validate' or 'validate'."
    }
    if (!params.ena_input) {
        error "--ena_input is required for the standalone ENA runner."
    }
    if (!params.ena_study?.toString()?.trim()) {
        error "--ena_study is required for Webin validation."
    }
    if (!params.outdir) {
        error "--outdir is required for the standalone ENA runner."
    }
    if (!(params.ena_validation_attempt?.toString() ==~ /[A-Za-z0-9._-]+/)) {
        error "--ena_validation_attempt may contain only letters, numbers, dots, underscores, and hyphens."
    }
    if (!secrets.WEBIN_USERNAME || !secrets.WEBIN_PASSWORD) {
        error "Nextflow secrets WEBIN_USERNAME and WEBIN_PASSWORD are required for Webin validation."
    }

    input_sheet = file(params.ena_input, checkIfExists: true)
    seen_keys = Collections.synchronizedSet(new HashSet<String>())
    ch_rows = channel.fromPath(input_sheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def rowLabel = "row for sample '${row.sample ?: 'unknown'}'"
            def sample = safeIdentifier(requiredValue(row, 'sample', rowLabel), 'sample', rowLabel)
            def prefix = safeIdentifier(requiredValue(row, 'mt_assembly_prefix', rowLabel), 'mt_assembly_prefix', rowLabel)
            def key = "${sample}\t${prefix}"
            if (!seen_keys.add(key)) {
                error "Duplicate sample/mt_assembly_prefix combination in ENA input: ${sample}/${prefix}"
            }
            tuple(row, [id: sample, mt_assembly_prefix: prefix])
        }
        .ifEmpty { error "ENA input CSV contains no data rows: ${input_sheet}" }

    if (params.ena_mode == 'convert_validate') {
        ch_prepared = ch_rows.map { row, meta ->
            def gbfValue = requiredValue(row, 'gbf', 0)
            def statusValue = requiredValue(row, 'table2asn_status', 0)
            if (!gbfValue.toLowerCase().endsWith('.gbf')) {
                error "GBF input for ${meta.id}/${meta.mt_assembly_prefix} must end in .gbf: ${gbfValue}"
            }
            def gbf = resolveInputPath(input_sheet, gbfValue, 'gbf', 0)
            def table2asnStatus = resolveInputPath(input_sheet, statusValue, 'table2asn_status', 0)
            def gate = table2asnGate(table2asnStatus)
            tuple(meta, gbf, table2asnStatus, gate.status, gate.reason)
        }

        ch_conversion_input = ch_prepared
            .filter { _meta, _gbf, _statusFile, gateStatus, _reason -> gateStatus == 'PASS' }
            .map { meta, gbf, _statusFile, _gateStatus, _reason -> tuple(meta, gbf) }

        ENA_FLATFILE(ch_conversion_input)
        WEBIN_VALIDATE(
            ENA_FLATFILE.out.embl_file,
            params.ena_study.toString().trim(),
            params.ena_validation_attempt.toString()
        )

        ch_validation_files = ch_prepared
            .map { meta, _gbf, statusFile, _gateStatus, _reason -> tuple(meta, statusFile) }
            .mix(ENA_FLATFILE.out.status)
            .mix(ENA_FLATFILE.out.checks)
            .mix(ENA_FLATFILE.out.embl_file)
            .mix(WEBIN_VALIDATE.out.status)
            .mix(WEBIN_VALIDATE.out.manifest)
        ch_multiqc_files = ch_prepared
            .map { _meta, _gbf, statusFile, _gateStatus, _reason -> statusFile }
            .mix(ENA_FLATFILE.out.status.map { _meta, statusFile -> statusFile })
            .mix(WEBIN_VALIDATE.out.status.map { _meta, statusFile -> statusFile })
    }
    else {
        ch_preflight_input = ch_rows.map { row, meta ->
            def emblValue = requiredValue(row, 'embl', 0)
            if (!emblValue.toLowerCase().endsWith('.embl.gz')) {
                error "EMBL input for ${meta.id}/${meta.mt_assembly_prefix} must end in .embl.gz: ${emblValue}"
            }
            def embl = resolveInputPath(input_sheet, emblValue, 'embl', 0)
            tuple(meta, embl)
        }
        ENA_EMBL_PREFLIGHT(ch_preflight_input)
        WEBIN_VALIDATE(
            ENA_EMBL_PREFLIGHT.out.embl_file,
            params.ena_study.toString().trim(),
            params.ena_validation_attempt.toString()
        )

        ch_validation_files = ENA_EMBL_PREFLIGHT.out.status
            .mix(ENA_EMBL_PREFLIGHT.out.checks)
            .mix(ENA_EMBL_PREFLIGHT.out.embl_file)
            .mix(WEBIN_VALIDATE.out.status)
            .mix(WEBIN_VALIDATE.out.manifest)
        ch_multiqc_files = ENA_EMBL_PREFLIGHT.out.status
            .map { _meta, statusFile -> statusFile }
            .mix(WEBIN_VALIDATE.out.status.map { _meta, statusFile -> statusFile })
    }

    ch_validation_inputs = ch_validation_files
        .map { meta, validationFile -> tuple(meta.mt_assembly_prefix, meta, validationFile) }
        .groupTuple(by: 0)
        .map { _prefix, metas, validationFiles -> tuple(metas[0], validationFiles.flatten()) }

    validation_settings = [
        ena_study: params.ena_study.toString().trim(),
        validation_mode: params.ena_mode,
        validation_attempt: params.ena_validation_attempt,
        webin_requested: true,
        workflow_run_name: workflow.runName ?: '',
        workflow_session_id: workflow.sessionId?.toString() ?: '',
        pipeline_revision: workflow.revision ?: workflow.commitId ?: ''
    ]
    ENA_VALIDATION_RESULT(ch_validation_inputs, validation_settings)
    ENA_VALIDATION_SUMMARY(ENA_VALIDATION_RESULT.out.record.map { _meta, record -> record }.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ENA_VALIDATION_SUMMARY.out.multiqc)

    if (!params.skip_upload_results && params.sql_config) {
        sql_config = file(params.sql_config, checkIfExists: true)
        UPLOAD_ENA_RESULTS(
            ENA_VALIDATION_RESULT.out.record,
            Channel.empty(),
            sql_config
        )
        ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_ENA_RESULTS.out.multiqc_files)
    }
    else if (!params.skip_upload_results && !params.sql_config) {
        log.warn 'ENA validation results were not uploaded: provide --sql_config or set --skip_upload_results true for file-only operation.'
    }

    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).toList()
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true).toList() : Channel.empty().toList()
    ch_multiqc_logo = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true).toList() : Channel.empty().toList()
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        [],
        []
    )

}
