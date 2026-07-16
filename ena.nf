#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENA_FLATFILE       } from './modules/local/genome_qc/ena_flatfile'
include { ENA_EMBL_PREFLIGHT } from './modules/local/genome_qc/ena_embl_preflight'
include { WEBIN_VALIDATE     } from './modules/local/genome_qc/webin_validate'

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

process ENA_RUN_SUMMARY {
    tag "${mode}:${validation_attempt}"
    label 'process_single'

    input:
    val rows
    val mode
    val validation_attempt

    output:
    path 'ena_run_summary.tsv', emit: summary

    script:
    def body = rows.sort().join('\n')
    """
    printf '%s\n' 'sample\tmt_assembly_prefix\tmode\ttable2asn_gate\tpreflight_status\tconversion_status\twebin_status\twebin_reason\tsubmission_ready\tvalidation_attempt' > ena_run_summary.tsv
    printf '%s\n' '${body}' >> ena_run_summary.tsv
    """
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
            tuple(meta, gbf, gate.status, gate.reason)
        }

        ch_gate_events = ch_prepared.map { meta, _gbf, gateStatus, reason ->
            tuple(meta.mt_assembly_prefix, meta.id, 'table2asn', gateStatus, reason)
        }
        ch_conversion_input = ch_prepared
            .filter { _meta, _gbf, gateStatus, _reason -> gateStatus == 'PASS' }
            .map { meta, gbf, _gateStatus, _reason -> tuple(meta, gbf) }

        ENA_FLATFILE(ch_conversion_input)
        WEBIN_VALIDATE(
            ENA_FLATFILE.out.embl_file,
            params.ena_study.toString().trim(),
            params.ena_validation_attempt.toString()
        )

        ch_conversion_events = ENA_FLATFILE.out.status.map { meta, statusFile ->
            def parsed = readNamedStatus(statusFile)
            tuple(meta.mt_assembly_prefix, meta.id, 'conversion', parsed.status, parsed.reason)
        }
        ch_webin_events = WEBIN_VALIDATE.out.status.map { meta, statusFile ->
            def parsed = readNamedStatus(statusFile)
            tuple(meta.mt_assembly_prefix, meta.id, 'webin', parsed.status, parsed.reason)
        }
        ch_events = ch_gate_events.mix(ch_conversion_events, ch_webin_events)
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
        ch_input_events = ch_preflight_input.map { meta, _embl ->
            tuple(meta.mt_assembly_prefix, meta.id, 'input', 'PASS', 'provided_embl')
        }

        ENA_EMBL_PREFLIGHT(ch_preflight_input)
        WEBIN_VALIDATE(
            ENA_EMBL_PREFLIGHT.out.embl_file,
            params.ena_study.toString().trim(),
            params.ena_validation_attempt.toString()
        )

        ch_preflight_events = ENA_EMBL_PREFLIGHT.out.status.map { meta, statusFile ->
            def parsed = readNamedStatus(statusFile)
            tuple(meta.mt_assembly_prefix, meta.id, 'preflight', parsed.status, parsed.reason)
        }
        ch_webin_events = WEBIN_VALIDATE.out.status.map { meta, statusFile ->
            def parsed = readNamedStatus(statusFile)
            tuple(meta.mt_assembly_prefix, meta.id, 'webin', parsed.status, parsed.reason)
        }
        ch_events = ch_input_events.mix(ch_preflight_events, ch_webin_events)
    }

    ch_summary_rows = ch_events
        .groupTuple(by: [0, 1])
        .map { prefix, sample, stages, statuses, reasons ->
            def events = [:]
            stages.eachWithIndex { stage, index ->
                events[stage] = [status: statuses[index], reason: reasons[index]]
            }
            def gate = params.ena_mode == 'convert_validate' ? (events.table2asn?.status ?: 'NOT_RUN') : 'NOT_APPLICABLE'
            def preflight = params.ena_mode == 'validate' ? (events.preflight?.status ?: 'NOT_RUN') : 'NOT_APPLICABLE'
            def conversion = params.ena_mode == 'convert_validate' ?
                (events.conversion?.status ?: (gate == 'PASS' ? 'NOT_RUN' : 'SKIPPED_TABLE2ASN')) :
                'NOT_APPLICABLE'
            def webin = events.webin?.status ?: 'NOT_RUN'
            def webinReason = events.webin?.reason ?: 'not_run'
            def ready = webin == 'PASS' ? 'true' : 'false'
            [sample, prefix, params.ena_mode, gate, preflight, conversion, webin, webinReason, ready, params.ena_validation_attempt].join('\t')
        }

    ENA_RUN_SUMMARY(ch_summary_rows.collect(), params.ena_mode, params.ena_validation_attempt)

}
