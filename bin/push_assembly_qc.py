#!/usr/bin/env python3
"""Store each assembly's QC verdict from mitogenome_assembly_summary_mqc.tsv.

MITOGENOME_ASSEMBLY_SUMMARY decides, per assembly, whether it is complete, needs manual
review, failed, or was superseded by a reseed or collapse, and why; whether it carries a
NUMT signal or a frameshift; and which reference MitoHiFi or OATK assembled against. Until
genomes_db_review migration 068 that verdict lived only in the run's MultiQC table.

Phase R: the verdict is a row of its own in mito_assembly_qc, one per assembly, keyed on the
assembly_id of the mito_assembly row push_mtdna_assm_results.py created. It never creates an
assembly: a verdict for one the database does not hold is reported as no_row and left alone. Every column is rewritten each time, so a verdict that
changes on a re-run replaces the old one, and qc_updated_at says when.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Optional, Union

try:
    import psycopg2
except ImportError:  # Allows unit tests to inject a connection factory.
    psycopg2 = None

# load_db_config is identical across every push script in bin/; reuse rather than copy.
from push_ena_validation_results import load_db_config


STATUSES = {"complete", "manual_review", "failed", "superseded"}

UPDATE_QUERY = """
    INSERT INTO mito_assembly_qc (assembly_id, qc_status, manual_review_reason, numt_flag,
                                  frameshift_flag, reference_species, reference_accession,
                                  qc_updated_at)
    SELECT a.assembly_id, %(qc_status)s, %(manual_review_reason)s, %(numt_flag)s,
           %(frameshift_flag)s, %(reference_species)s, %(reference_accession)s, now()
      FROM mito_assembly a
     WHERE a.og_id = %(og_id)s
       AND a.tech = %(tech)s
       AND a.seq_date = %(seq_date)s
       AND a.code = %(code)s
    ON CONFLICT (assembly_id) DO UPDATE
       SET qc_status = EXCLUDED.qc_status,
           manual_review_reason = EXCLUDED.manual_review_reason,
           numt_flag = EXCLUDED.numt_flag,
           frameshift_flag = EXCLUDED.frameshift_flag,
           reference_species = EXCLUDED.reference_species,
           reference_accession = EXCLUDED.reference_accession,
           qc_updated_at = now()
"""


def _flag(value: Optional[str]) -> Optional[bool]:
    """mitogenome_assembly_summary.py writes "true", "false" or blank (not determined)."""
    text = (value or "").strip().lower()
    if text == "true":
        return True
    if text == "false":
        return False
    if text == "":
        return None
    raise ValueError(f"unexpected flag value {value!r}")


def _text(value: Optional[str]) -> Optional[str]:
    text = (value or "").strip()
    return text or None


def read_verdicts(path: Union[str, Path]) -> tuple[list[dict], list[str]]:
    """Parse the summary table into one record per keyable assembly.

    Returns (records, skipped). An assembly prefix without the four dot fields
    og_id.tech.seq_date.code cannot be matched to a row and is skipped, by name.
    """
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    records, skipped = [], []
    for row in rows:
        prefix = (row.get("assembly_prefix") or "").strip()
        parts = prefix.split(".")
        if len(parts) != 4 or not all(parts):
            skipped.append(prefix or "<blank assembly_prefix>")
            continue
        status = (row.get("status") or "").strip()
        if status not in STATUSES:
            raise ValueError(f"{prefix}: unexpected status {status!r}")
        og_id, tech, seq_date, code = parts
        records.append({
            "og_id": og_id, "tech": tech, "seq_date": seq_date, "code": code,
            "qc_status": status,
            "manual_review_reason": _text(row.get("manual_review_reason")),
            "numt_flag": _flag(row.get("numt_flag")),
            "frameshift_flag": _flag(row.get("frameshift_flag")),
            "reference_species": _text(row.get("reference_species")),
            "reference_accession": _text(row.get("reference_accession")),
        })
    return records, skipped


def push_verdicts(records: list[dict], db_config: dict, dry_run: bool = False,
                  connect=None) -> dict[str, list[str]]:
    """Write every record in one transaction; return prefixes by outcome."""
    connect = connect or (psycopg2.connect if psycopg2 is not None else None)
    if connect is None:
        raise RuntimeError("psycopg2 is required for PostgreSQL upload")
    outcome: dict[str, list[str]] = {"set": [], "no_row": []}
    connection = connect(**db_config)
    try:
        with connection.cursor() as cursor:
            for record in records:
                cursor.execute(UPDATE_QUERY, record)
                prefix = ".".join(record[k] for k in ("og_id", "tech", "seq_date", "code"))
                outcome["set" if cursor.rowcount == 1 else "no_row"].append(prefix)
        if dry_run:
            connection.rollback()
        else:
            connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()
    return outcome


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("config_file")
    parser.add_argument("summary_tsv", help="mitogenome_assembly_summary_mqc.tsv")
    parser.add_argument("--dry-run", action="store_true",
                        help="run every update, report what would change, then roll back")
    args = parser.parse_args()

    try:
        records, skipped = read_verdicts(args.summary_tsv)
        outcome = push_verdicts(records, load_db_config(args.config_file), dry_run=args.dry_run)
    except Exception as exc:
        print(f"❌ Error: {exc}", file=sys.stderr)
        return 1

    verb = "would set" if args.dry_run else "set"
    print(f"✅ {verb} the QC verdict on {len(outcome['set'])} assembly row(s)")
    for prefix in outcome["no_row"]:
        print(f"⚠️  no mitogenome_data row for {prefix}; verdict not stored")
    for prefix in skipped:
        print(f"⚠️  skipped {prefix}: not og_id.tech.seq_date.code")
    return 0


if __name__ == "__main__":
    sys.exit(main())
