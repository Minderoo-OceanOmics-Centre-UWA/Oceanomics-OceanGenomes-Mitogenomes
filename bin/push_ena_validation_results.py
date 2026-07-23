#!/usr/bin/env python3
"""Insert a normalized ENA validation result into PostgreSQL."""

from __future__ import annotations

import argparse
import configparser
import csv
import sys
from pathlib import Path
from typing import Union

try:
    import psycopg2
except ImportError:  # Allows unit tests to inject a connection factory.
    psycopg2 = None


INTEGER_COLUMNS = {
    "reject_count", "error_count", "warning_count", "info_count",
    "fatal_discrepancy_count", "nostop_count", "conversion_exit",
    "preflight_exit", "webin_exit", "flatfile_size", "manifest_size",
}

INSERT_COLUMNS = [
    "assembly_prefix", "og_id", "tech", "seq_date", "code", "ena_study",
    "validation_mode", "validation_attempt", "table2asn_status",
    "reject_count", "error_count", "warning_count", "info_count",
    "fatal_discrepancy_count", "nostop_count", "blocking_codes", "warning_codes",
    "conversion_status", "conversion_reason", "conversion_exit",
    "preflight_status", "preflight_reason", "preflight_exit", "webin_status",
    "webin_reason", "webin_exit", "submission_ready", "flatfile_name",
    "flatfile_sha256", "flatfile_size", "manifest_name", "manifest_sha256",
    "manifest_size", "workflow_run_name", "workflow_session_id",
    "pipeline_revision", "result_digest",
]


def load_db_config(path: Union[str, Path]) -> dict[str, object]:
    parser = configparser.ConfigParser()
    if not parser.read(path) or not parser.has_section("postgres"):
        raise ValueError(f"Missing [postgres] configuration in {path}")
    required = ["dbname", "user", "password", "host", "port"]
    missing = [key for key in required if not parser.has_option("postgres", key)]
    if missing:
        raise ValueError(f"Missing PostgreSQL settings: {', '.join(missing)}")
    return {
        "dbname": parser.get("postgres", "dbname"),
        "user": parser.get("postgres", "user"),
        "password": parser.get("postgres", "password"),
        "host": parser.get("postgres", "host"),
        "port": parser.getint("postgres", "port"),
    }


def read_record(path: Union[str, Path]) -> dict[str, object]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one ENA validation row in {path}")
    missing = [column for column in INSERT_COLUMNS if column not in rows[0]]
    if missing:
        raise ValueError(f"Missing ENA validation columns: {', '.join(missing)}")
    record: dict[str, object] = {}
    for column in INSERT_COLUMNS:
        value = (rows[0].get(column) or "").strip()
        if column in INTEGER_COLUMNS:
            record[column] = int(value) if value else None
        elif column == "submission_ready":
            record[column] = value.lower() in {"true", "1", "yes"}
        elif column == "ena_study":
            record[column] = value
        else:
            record[column] = value or None
    for required in ("assembly_prefix", "og_id", "validation_mode", "validation_attempt",
                     "table2asn_status", "conversion_status", "preflight_status",
                     "webin_status", "result_digest"):
        if record.get(required) is None:
            raise ValueError(f"Required ENA validation value is empty: {required}")
    return record


KEY_COLUMNS = {"assembly_prefix", "ena_study", "validation_attempt"}


def upload_record(record: dict[str, object], db_config: dict[str, object], connect=None) -> str:
    """Insert or overwrite the row for this (assembly_prefix, ena_study, validation_attempt).

    A rerun under the same key overwrites the previous attempt so failed
    attempts don't pile up history rows. Once a row's submission_ready is
    true it is frozen: the WHERE clause on the DO UPDATE skips the overwrite,
    so a later failing rerun can never clobber a recorded success.

    Returns "inserted" (first row for this key), "updated" (overwrote a
    not-yet-ready attempt), or "locked" (a submission-ready row already
    exists and was left untouched).
    """
    connect = connect or (psycopg2.connect if psycopg2 is not None else None)
    if connect is None:
        raise RuntimeError("psycopg2 is required for PostgreSQL upload")
    update_columns = [column for column in INSERT_COLUMNS if column not in KEY_COLUMNS]
    set_clause = ", ".join(f"{column} = EXCLUDED.{column}" for column in update_columns)
    placeholders = ", ".join(f"%({column})s" for column in INSERT_COLUMNS)
    columns = ", ".join(INSERT_COLUMNS)
    query = f"""
        INSERT INTO ena_validation_attempts ({columns})
        VALUES ({placeholders})
        ON CONFLICT (assembly_prefix, ena_study, validation_attempt)
        DO UPDATE SET
            {set_clause},
            recorded_at = CURRENT_TIMESTAMP,
            attempt_count = ena_validation_attempts.attempt_count + 1
        WHERE ena_validation_attempts.submission_ready = FALSE
        RETURNING id, (xmax = 0) AS inserted
    """
    connection = connect(**db_config)
    try:
        with connection.cursor() as cursor:
            cursor.execute(query, record)
            row = cursor.fetchone()
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()
    if row is None:
        return "locked"
    return "inserted" if row[1] else "updated"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file")
    parser.add_argument("validation_record")
    args = parser.parse_args()
    try:
        record = read_record(args.validation_record)
        result = upload_record(record, load_db_config(args.config_file))
        if result == "inserted":
            print(f"✅ Success: inserted ENA validation attempt for {record['assembly_prefix']}")
        elif result == "updated":
            print(
                f"🔁 Updated ENA validation attempt for {record['assembly_prefix']} "
                "(previous attempt was not submission-ready)"
            )
        else:
            print(
                f"🔒 Submission-ready record already exists for {record['assembly_prefix']}; "
                "new attempt not recorded"
            )
        return 0
    except Exception as error:
        print(f"❌ Database error: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
