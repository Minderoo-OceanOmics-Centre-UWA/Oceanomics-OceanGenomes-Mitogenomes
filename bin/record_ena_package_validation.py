#!/usr/bin/env python3
"""Persist a genome-context Webin validation result without submitting it."""

from __future__ import annotations

import argparse
import configparser
import csv
import sys
from pathlib import Path

try:
    import psycopg2
except ImportError:
    psycopg2 = None


def load_config(path: Path) -> dict[str, object]:
    parser = configparser.ConfigParser()
    parser.read(path)
    if not parser.has_section("postgres"):
        raise ValueError(f"Missing [postgres] section in {path}")
    return {
        "dbname": parser.get("postgres", "dbname"),
        "user": parser.get("postgres", "user"),
        "password": parser.get("postgres", "password"),
        "host": parser.get("postgres", "host"),
        "port": parser.getint("postgres", "port"),
    }


def read_status(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one validation row in {path}")
    required = {"full_seqid", "service", "status", "reason"}
    missing = required - rows[0].keys()
    if missing:
        raise ValueError(f"Missing validation columns: {', '.join(sorted(missing))}")
    return {key: (value or "").strip() for key, value in rows[0].items()}


def persist(connection, row: dict[str, str]) -> None:
    service = row["service"]
    if service not in {"test", "production"}:
        raise ValueError(f"Invalid Webin service: {service}")
    column = f"webin_{service}_status"
    with connection.cursor() as cursor:
        cursor.execute(
            f"""
            UPDATE ena_candidate_packages
            SET {column} = %s, updated_at = CURRENT_TIMESTAMP
            WHERE full_seqid = %s
            RETURNING assembly_prefix, og_id, ena_study_accession
            """,
            (row["status"], row["full_seqid"]),
        )
        package = cursor.fetchone()
        if package is None:
            raise ValueError(
                f"Candidate package is not registered: {row['full_seqid']}"
            )
        assembly_prefix, og_id, study = package
        if service == "production":
            cursor.execute(
                """
                UPDATE ena_validation_attempts validation
                SET webin_production_status = %s,
                    webin_reason = %s,
                    webin_exit = %s,
                    submission_ready = (
                        %s = 'PASS'
                        AND EXISTS (
                            SELECT 1
                            FROM ena_submission_selections selection
                            WHERE selection.ena_study_accession = %s
                              AND selection.og_id = %s
                              AND selection.selected_full_seqid = %s
                              AND selection.selection_status = 'SELECTED'
                        )
                    ),
                    overall_status = CASE
                        WHEN %s = 'PASS' THEN 'PRODUCTION_VALIDATED'
                        ELSE %s
                    END,
                    recorded_at = CURRENT_TIMESTAMP
                WHERE validation.assembly_prefix = %s
                  AND validation.ena_study = %s
                """,
                (
                    row["status"],
                    row["reason"],
                    int(row["webin_exit"]) if row.get("webin_exit") else None,
                    row["status"],
                    study,
                    og_id,
                    row["full_seqid"],
                    row["status"],
                    row["status"],
                    assembly_prefix,
                    study,
                ),
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--status", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        if psycopg2 is None:
            raise RuntimeError("psycopg2 is required")
        row = read_status(Path(args.status))
        with psycopg2.connect(**load_config(Path(args.config))) as connection:
            persist(connection, row)
        Path(args.output).write_text(
            "full_seqid\tservice\tstatus\trecorded\n"
            f"{row['full_seqid']}\t{row['service']}\t{row['status']}\ttrue\n"
        )
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
