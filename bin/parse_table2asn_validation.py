#!/usr/bin/env python3
"""Normalise table2asn validator/discrepancy output without failing a batch."""

import argparse
import csv
import re
from pathlib import Path


VALIDATOR_RE = re.compile(r"^(REJECT|ERROR|WARNING|INFO):\s*[^[]*\[([^]]+)]\s*(.*)$")
FATAL_RE = re.compile(r"(?:^|\s)FATAL(?::|\s|$)")
FATAL_CODE_RE = re.compile(r"FATAL[:\s]+([^:\s]+)")


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
        else:
            severity, code, message = "INFO", "UNPARSED", raw
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
            findings.append((path.name, "FATAL", code, raw.replace("\t", " ")))
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
