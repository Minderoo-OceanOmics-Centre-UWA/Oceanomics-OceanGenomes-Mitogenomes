#!/usr/bin/env python3
"""Normalize ENA submission-gate outputs and compile batch summaries."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
from pathlib import Path


RECORD_COLUMNS = [
    "assembly_prefix", "og_id", "tech", "seq_date", "code",
    "ena_study", "validation_mode", "validation_attempt",
    "table2asn_status", "reject_count", "error_count", "warning_count",
    "info_count", "fatal_discrepancy_count", "nostop_count",
    "blocking_codes", "warning_codes",
    "conversion_status", "conversion_reason", "conversion_exit",
    "preflight_status", "preflight_reason", "preflight_exit",
    "webin_status", "webin_reason", "webin_exit", "submission_ready",
    "package_status", "local_package_status", "webin_test_status",
    "webin_production_status", "webin_error_count", "webin_warning_count",
    "webin_error_codes", "webin_cli_version", "webin_test_report_path",
    "webin_production_report_path", "overall_status", "package_digest",
    "flatfile_name", "flatfile_sha256", "flatfile_size",
    "manifest_name", "manifest_sha256", "manifest_size",
    "workflow_run_name", "workflow_session_id", "pipeline_revision",
    "result_digest",
]

MQC_COLUMNS = [
    "assembly_prefix", "og_id", "validation_mode", "ena_study",
    "table2asn_status", "conversion_status", "preflight_status",
    "webin_status", "webin_reason", "submission_ready",
    "package_status", "local_package_status", "webin_test_status",
    "webin_production_status", "overall_status",
    "validation_attempt",
]

DIGEST_COLUMNS = [
    column for column in RECORD_COLUMNS
    if column not in {"result_digest", "workflow_run_name", "workflow_session_id", "pipeline_revision"}
]


def read_tsv_row(path: Path) -> dict[str, str]:
    try:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return {}
    if not rows:
        return {}
    return {str(key): "" if value is None else str(value).strip() for key, value in rows[0].items()}


def read_status(path: Path) -> dict[str, str]:
    row = read_tsv_row(path)
    if not row.get("status"):
        return {"status": "MALFORMED_STATUS", "reason": "missing_status_row"}
    return row


def file_metadata(path: Path | None) -> tuple[str, str, str]:
    if path is None or not path.is_file():
        return "", "", ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return path.name, digest.hexdigest(), str(path.stat().st_size)


def identity_from_prefix(prefix: str, og_id: str) -> tuple[str, str, str, str]:
    parts = prefix.split(".", 3)
    inferred_og = parts[0] if parts and parts[0] else og_id
    tech = parts[1] if len(parts) > 1 else ""
    seq_date = parts[2] if len(parts) > 2 else ""
    code = parts[3] if len(parts) > 3 else ""
    return inferred_og or og_id, tech, seq_date, code


def digest_record(record: dict[str, str]) -> str:
    payload = "\n".join(f"{column}={record.get(column, '')}" for column in DIGEST_COLUMNS)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _first(paths: list[Path], suffix: str) -> Path | None:
    return next((path for path in paths if path.name.endswith(suffix)), None)


def build_record(
    input_paths: list[Path], *, assembly_prefix: str, og_id: str,
    ena_study: str, validation_mode: str, validation_attempt: str,
    webin_requested: bool, workflow_run_name: str = "",
    workflow_session_id: str = "", pipeline_revision: str = "",
) -> dict[str, str]:
    table_path = _first(input_paths, ".table2asn_status.tsv")
    conversion_path = _first(input_paths, ".ena_conversion_status.tsv")
    preflight_path = _first(input_paths, ".ena_preflight_status.tsv")
    preflight_check_path = _first(input_paths, ".ena_preflight_check.tsv")
    webin_path = _first(input_paths, ".webin_status.tsv")
    webin_test_path = _first(input_paths, ".webin_test_status.tsv")
    webin_production_path = _first(input_paths, ".webin_production_status.tsv")
    local_package_path = _first(input_paths, ".local_validation.tsv")
    package_metadata_path = _first(input_paths, ".package_metadata.json")
    flatfile_path = next((p for p in input_paths if p.name.endswith(".embl.gz")), None)
    manifest_path = _first(input_paths, ".webin_manifest.txt")

    table = read_status(table_path) if table_path else {}
    conversion = read_status(conversion_path) if conversion_path else {}
    preflight = read_status(preflight_path) if preflight_path else {}
    preflight_check = read_tsv_row(preflight_check_path) if preflight_check_path else {}
    webin = read_status(webin_path) if webin_path else {}
    webin_test = read_status(webin_test_path) if webin_test_path else {}
    webin_production = read_status(webin_production_path) if webin_production_path else {}
    local_package = read_status(local_package_path) if local_package_path else {}
    package_metadata = {}
    if package_metadata_path:
        try:
            package_metadata = json.loads(package_metadata_path.read_text())
        except (OSError, json.JSONDecodeError):
            package_metadata = {}

    if validation_mode == "validate":
        table_status = "NOT_APPLICABLE"
        conversion_status = "NOT_APPLICABLE"
        preflight_status = preflight.get("status", "NOT_RUN")
        webin_status = webin.get("status", "NOT_RUN")
    else:
        table_status = table.get("status", "MALFORMED_STATUS")
        conversion_status = conversion.get(
            "status", "NOT_RUN" if table_status == "PASS" else "SKIPPED_TABLE2ASN"
        )
        preflight_status = "NOT_APPLICABLE"
        if not webin_requested:
            webin_status = "NOT_REQUESTED"
        else:
            webin_status = webin.get("status", "NOT_RUN")

    if validation_mode == "validate" and not webin_requested:
        webin_status = "NOT_REQUESTED"

    if conversion:
        conversion_reason = conversion.get("reason", "unspecified")
    elif validation_mode == "validate":
        conversion_reason = "not_applicable"
    elif table_status != "PASS":
        conversion_reason = "table2asn_blocked"
    else:
        conversion_reason = "not_run"

    if webin:
        webin_reason = webin.get("reason", "unspecified")
    elif not webin_requested:
        webin_reason = "not_requested"
    elif ((validation_mode == "validate" and preflight_status != "PASS") or
          (validation_mode != "validate" and conversion_status != "PASS")):
        webin_reason = "upstream_gate_not_passed"
    else:
        webin_reason = "not_run"

    inferred_og, tech, seq_date, code = identity_from_prefix(assembly_prefix, og_id)
    flat_name, flat_hash, flat_size = file_metadata(flatfile_path)
    manifest_name, manifest_hash, manifest_size = file_metadata(manifest_path)
    gates_pass = (
        preflight_status == "PASS" if validation_mode == "validate"
        else table_status == "PASS" and conversion_status == "PASS"
    )
    package_status = str(package_metadata.get("package_status") or "NOT_RUN")
    local_package_status = local_package.get(
        "status", str(package_metadata.get("local_validation_status") or "NOT_RUN")
    )
    webin_test_status = webin_test.get("status", "NOT_RUN")
    webin_production_status = webin_production.get("status", "NOT_RUN")
    # Selection for submission is deliberately decided outside this per-candidate
    # collator, so a candidate can never become submission-ready here.
    ready = "false"
    if webin_production_status == "PASS":
        overall_status = "PRODUCTION_VALIDATED"
    elif webin_production_status not in {"NOT_RUN", "NOT_REQUESTED"}:
        overall_status = webin_production_status
    elif webin_test_status == "PASS":
        overall_status = "WEBIN_TEST_PASS"
    elif webin_test_status not in {"NOT_RUN", "NOT_REQUESTED"}:
        overall_status = webin_test_status
    elif package_status == "READY" and local_package_status == "PASS":
        overall_status = "LOCAL_PACKAGE_READY"
    elif package_status != "NOT_RUN":
        overall_status = package_status
    elif gates_pass and webin_status == "PASS":
        overall_status = "LEGACY_SEQUENCE_VALIDATED"
    else:
        overall_status = "NOT_RUN"

    record = {column: "" for column in RECORD_COLUMNS}
    record.update({
        "assembly_prefix": assembly_prefix,
        "og_id": inferred_og,
        "tech": tech,
        "seq_date": seq_date,
        "code": code,
        "ena_study": ena_study,
        "validation_mode": validation_mode,
        "validation_attempt": validation_attempt,
        "table2asn_status": table_status,
        "reject_count": table.get("reject_count", ""),
        "error_count": table.get("error_count", ""),
        "warning_count": table.get("warning_count", ""),
        "info_count": table.get("info_count", ""),
        "fatal_discrepancy_count": table.get("fatal_discrepancy_count", ""),
        "nostop_count": table.get("nostop_count", ""),
        "blocking_codes": table.get("blocking_codes", ""),
        "warning_codes": table.get("warning_codes", ""),
        "conversion_status": conversion_status,
        "conversion_reason": conversion_reason,
        "conversion_exit": conversion.get("seqret_exit", ""),
        "preflight_status": preflight_status,
        "preflight_reason": preflight.get("reason", "not_applicable" if validation_mode != "validate" else "not_run"),
        "preflight_exit": preflight_check.get("gzip_exit", ""),
        "webin_status": webin_status,
        "webin_reason": webin_reason,
        "webin_exit": webin.get("webin_exit", ""),
        "submission_ready": ready,
        "package_status": package_status,
        "local_package_status": local_package_status,
        "webin_test_status": webin_test_status,
        "webin_production_status": webin_production_status,
        "webin_error_count": webin_production.get(
            "error_count", webin_test.get("error_count", "")
        ),
        "webin_warning_count": webin_production.get(
            "warning_count", webin_test.get("warning_count", "")
        ),
        "webin_error_codes": webin_production.get(
            "error_codes", webin_test.get("error_codes", "")
        ),
        "webin_cli_version": webin_production.get(
            "webin_cli_version", webin_test.get("webin_cli_version", "")
        ),
        "webin_test_report_path": webin_test.get("report_path", ""),
        "webin_production_report_path": webin_production.get("report_path", ""),
        "overall_status": overall_status,
        "package_digest": str(package_metadata.get("package_digest") or ""),
        "flatfile_name": flat_name,
        "flatfile_sha256": flat_hash,
        "flatfile_size": flat_size,
        "manifest_name": manifest_name,
        "manifest_sha256": manifest_hash,
        "manifest_size": manifest_size,
        "workflow_run_name": workflow_run_name,
        "workflow_session_id": workflow_session_id,
        "pipeline_revision": pipeline_revision,
    })
    record["result_digest"] = digest_record(record)
    return record


def write_rows(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_records(paths: list[Path]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="") as handle:
            records.extend(csv.DictReader(handle, delimiter="\t"))
    return sorted(records, key=lambda row: (row.get("assembly_prefix", ""), row.get("validation_attempt", "")))


def expand_patterns(patterns: list[str]) -> list[Path]:
    return [Path(match) for pattern in patterns for match in glob.glob(pattern)]


def command_record(args: argparse.Namespace) -> None:
    paths = expand_patterns(args.input)
    record = build_record(
        paths,
        assembly_prefix=args.assembly_prefix,
        og_id=args.og_id,
        ena_study=args.ena_study,
        validation_mode=args.validation_mode,
        validation_attempt=args.validation_attempt,
        webin_requested=args.webin_requested,
        workflow_run_name=args.workflow_run_name,
        workflow_session_id=args.workflow_session_id,
        pipeline_revision=args.pipeline_revision,
    )
    write_rows(Path(args.output), RECORD_COLUMNS, [record])


def command_summary(args: argparse.Namespace) -> None:
    paths = expand_patterns(args.input)
    rows = read_records(paths)
    write_rows(Path(args.output), MQC_COLUMNS, rows)
    if args.run_summary:
        write_rows(Path(args.run_summary), MQC_COLUMNS, rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    record = subparsers.add_parser("record")
    record.add_argument("--input", action="append", required=True)
    record.add_argument("--output", required=True)
    record.add_argument("--assembly-prefix", required=True)
    record.add_argument("--og-id", required=True)
    record.add_argument("--ena-study", default="")
    record.add_argument("--validation-mode", required=True, choices=["pipeline", "convert_validate", "validate"])
    record.add_argument("--validation-attempt", required=True)
    record.add_argument("--webin-requested", action="store_true")
    record.add_argument("--workflow-run-name", default="")
    record.add_argument("--workflow-session-id", default="")
    record.add_argument("--pipeline-revision", default="")
    record.set_defaults(func=command_record)

    summary = subparsers.add_parser("summary")
    summary.add_argument("--input", action="append", required=True)
    summary.add_argument("--output", required=True)
    summary.add_argument("--run-summary")
    summary.set_defaults(func=command_summary)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
