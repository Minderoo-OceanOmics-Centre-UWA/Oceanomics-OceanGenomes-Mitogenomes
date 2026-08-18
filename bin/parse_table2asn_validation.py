#!/usr/bin/env python3
"""Normalise table2asn validator/discrepancy output without failing a batch."""

import argparse
import csv
import re
from pathlib import Path


# table2asn writes the severity in mixed case ("Error: valid [SEQ_FEAT.StartCodon] ..."),
# so this must be case-insensitive. Matching only upper case silently routed every
# real finding into the UNPARSED/INFO fallback below, which made the gate report
# PASS for every sample in a run regardless of what the .val actually said.
VALIDATOR_RE = re.compile(r"^(REJECT|ERROR|WARNING|INFO):\s*[^[]*\[([^]]+)]\s*(.*)$",
                          re.IGNORECASE)
FATAL_RE = re.compile(r"(?:^|\s)FATAL(?::|\s|$)")
FATAL_CODE_RE = re.compile(r"FATAL[:\s]+([^:\s]+)")

# Discrepancy-report FATAL codes that are expected for organelle-only
# submissions and shouldn't quarantine a sample. NO_LOCUS_TAGS fires on EVERY
# record this pipeline produces, by design: locus tags are assigned by the
# downstream submission pipeline, so nothing here writes a /locus_tag and
# table2asn always sees a tag-free feature table. GenBank/ENA do not require
# locus tags on organelle genomes either. Do not promote this back to fatal --
# it would quarantine every sample in the run.
# MISSING_PROTEIN_ID fires on every sample because we deliberately omit
# protein_id (EMMA's placeholder UUID isn't a real accession; ENA/GenBank
# assign the real one at accessioning time), so it's expected, not an error.
ADVISORY_DISCREPANCY_CODES = {"NO_LOCUS_TAGS", "MISSING_PROTEIN_ID"}

# Validator codes that table2asn raises against NCBI submission rules which do
# not apply to this pipeline's ENA route. They are demoted to WARNING so they
# stay visible in the findings table without quarantining an otherwise good
# assembly. Everything not listed here keeps the severity table2asn assigned.
#
# LatLonWater/LatLonGeoLocName: NCBI cross-checks the coordinate against a land
#   polygon and complains that an offshore sample is "in water". ENA performs no
#   such check, and these fire on most OceanOmics samples by their nature.
# GeneXrefWithoutLocus: an artefact of deliberately shipping a locus_tag-free
#   feature table, for the same reason NO_LOCUS_TAGS is advisory above.
# OrganismIsUndefinedSpecies: NCBI wants a specific identifier appended to a
#   'Genus sp.' name. ENA requires the opposite -- webin-cli rejects anything
#   other than the bare submittable taxon (see species_name_utils
#   .normalise_open_nomenclature), and validates these samples today. Acting on
#   the NCBI advice here would break the ENA path.
ADVISORY_VALIDATOR_CODES = {
    "SEQ_DESCR.LatLonWater",
    "SEQ_DESCR.LatLonGeoLocName",
    "SEQ_FEAT.GeneXrefWithoutLocus",
    "SEQ_DESCR.OrganismIsUndefinedSpecies",
}


def unique_join(values):
    return ",".join(dict.fromkeys(values))


def parse_validator(path):
    findings = []
    if not path or not path.exists():
        return findings
    for raw in path.read_text(errors="replace").splitlines():
        if not raw.strip():
            continue
        match = VALIDATOR_RE.match(raw)
        if match:
            severity, code, message = match.groups()
            severity = severity.upper()
            if code in ADVISORY_VALIDATOR_CODES:
                severity = "WARNING"
        else:
            # An unrecognised line means table2asn changed its output format, not
            # that the record is clean. Surface it as a WARNING so it is visible
            # in the findings table rather than being buried among INFO rows.
            severity, code, message = "WARNING", "UNPARSED", raw
        findings.append((path.name, severity, code, message.replace("\t", " ")))
    return findings


def parse_discrepancy(path):
    findings = []
    if not path or not path.exists():
        return findings
    for raw in path.read_text(errors="replace").splitlines():
        if FATAL_RE.search(raw):
            code_match = FATAL_CODE_RE.search(raw)
            code = code_match.group(1) if code_match else "DISCREPANCY_FATAL"
            severity = "WARNING" if code in ADVISORY_DISCREPANCY_CODES else "FATAL"
            findings.append((path.name, severity, code, raw.replace("\t", " ")))
    return findings


def write_tsv(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--circular", required=True)
    parser.add_argument("--val", type=Path)
    parser.add_argument("--dr", type=Path)
    parser.add_argument("--findings", type=Path, required=True)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--qc-flags", type=Path, required=True)
    args = parser.parse_args()

    findings = parse_validator(args.val) + parse_discrepancy(args.dr)
    counts = {severity: 0 for severity in ("REJECT", "ERROR", "WARNING", "INFO", "FATAL")}
    for _source, severity, _code, _message in findings:
        counts[severity] = counts.get(severity, 0) + 1

    blocking_codes = [code for _, severity, code, _ in findings if severity in {"REJECT", "ERROR", "FATAL"}]
    warning_codes = [code for _, severity, code, _ in findings if severity == "WARNING"]
    nostop_messages = [message for _, _, code, message in findings if code == "SEQ_FEAT.NoStop"]
    nostop_features = []
    for message in nostop_messages:
        match = re.search(r"CDS:\s*(.*?)(?:\s+[<\[].*)?$", message)
        nostop_features.append(match.group(1).strip() if match else message.strip())

    gate_status = "FAIL_TABLE2ASN" if counts["REJECT"] or counts["ERROR"] or counts["FATAL"] else "PASS"
    write_tsv(
        args.findings,
        ["sample", "source", "severity", "code", "message"],
        [(args.sample, *finding) for finding in findings],
    )
    write_tsv(
        args.status,
        ["sample", "status", "reject_count", "error_count", "warning_count", "info_count", "fatal_discrepancy_count", "nostop_count", "blocking_codes", "warning_codes"],
        [[args.sample, gate_status, counts["REJECT"], counts["ERROR"], counts["WARNING"], counts["INFO"], counts["FATAL"], len(nostop_messages), unique_join(blocking_codes), unique_join(warning_codes)]],
    )
    write_tsv(
        args.qc_flags,
        ["sample", "circular", "status", "reject_count", "error_count", "warning_count", "fatal_discrepancy_count", "nostop_count", "nostop_features"],
        [[args.sample, args.circular, gate_status, counts["REJECT"], counts["ERROR"], counts["WARNING"], counts["FATAL"], len(nostop_messages), ";".join(nostop_features)]],
    )


if __name__ == "__main__":
    main()
