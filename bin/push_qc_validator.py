#!/usr/bin/env python3
"""Record the pipeline as the second species-ID validator in ``lca_validation``.

A mitogenome needs two validators signed off in ``lca_validation`` before it is
OK to submit.  ``species_validation.py`` fills the first slot (``validator =
'nf-core'``) when the nominal species ID is found in the BLAST results; the
second slot was filled by hand.

This script fills the second slot automatically for any sample that cleared
every QC gate in the pipeline -- table2asn, the EMBL flat-file conversion and
the ``ena-webin-cli -validate`` format check, i.e. ``submission_ready = true``
in the per-sample ENA validation record.  It reads the same
``<full_seqid>.ena_validation_result.tsv`` that ``push_ena_validation_results.py``
consumes, so it needs no input of its own.

The write NEVER overwrites an existing ``validator_2``.  A human reviewer who
already signed off stays recorded, and there is deliberately no ``--force``
escape hatch: the point of the column is that a second, independent validator
looked at the sample, so clobbering one is never the right move.  The row is
also never INSERTed -- with no ``lca_validation`` row there is no first
validator either, and a lone ``validator_2`` would be meaningless.
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

# load_db_config is identical across every push script in bin/; reuse rather
# than copy. bin/ is on PATH at runtime and sys.path[0] is this file's
# directory, so the sibling import resolves the same way species_validation.py
# imports species_name_utils.
from push_ena_validation_results import load_db_config


VALIDATOR_2 = "QCd-nf-core"

# The lca_validation primary key, all of which the ENA validation record
# already carries.
KEY_COLUMNS = ["og_id", "tech", "seq_date", "code", "annotation"]
IDENTITY_COLUMNS = ["full_seqid"] + KEY_COLUMNS

UPDATE_QUERY = """
    UPDATE lca_validation
       SET validator_2 = %(validator_2)s
     WHERE og_id = %(og_id)s
       AND tech = %(tech)s
       AND seq_date = %(seq_date)s
       AND code = %(code)s
       AND annotation = %(annotation)s
       AND (validator_2 IS NULL OR btrim(validator_2) = '')
    RETURNING validator_2
"""

EXISTING_QUERY = """
    SELECT validator, validator_2
      FROM lca_validation
     WHERE og_id = %(og_id)s
       AND tech = %(tech)s
       AND seq_date = %(seq_date)s
       AND code = %(code)s
       AND annotation = %(annotation)s
"""


def read_record(path: Union[str, Path]) -> dict:
    """Read the identity fields and the submission_ready verdict from the record."""
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one ENA validation row in {path}")
    missing = [column for column in IDENTITY_COLUMNS + ["submission_ready"]
               if column not in rows[0]]
    if missing:
        raise ValueError(f"Missing ENA validation columns: {', '.join(missing)}")
    record: dict = {
        column: (rows[0].get(column) or "").strip() for column in IDENTITY_COLUMNS
    }
    empty = [column for column in IDENTITY_COLUMNS if not record[column]]
    if empty:
        raise ValueError(f"Required ENA validation value is empty: {', '.join(empty)}")
    ready = (rows[0].get("submission_ready") or "").strip().lower()
    record["submission_ready"] = ready in {"true", "1", "yes"}
    return record


def set_qc_validator(record: dict, db_config: dict, connect=None) -> tuple[str, Optional[str]]:
    """Fill validator_2 for this sample if and only if the column is still blank.

    Returns ``(status, detail)`` where status is one of:

      * ``"set"``       -- the column was blank and now holds VALIDATOR_2
      * ``"preserved"`` -- the row already carried a validator_2 (``detail``)
      * ``"no_row"``    -- no lca_validation row exists for this key
    """
    connect = connect or (psycopg2.connect if psycopg2 is not None else None)
    if connect is None:
        raise RuntimeError("psycopg2 is required for PostgreSQL upload")
    params = {column: record[column] for column in KEY_COLUMNS}
    params["validator_2"] = VALIDATOR_2
    connection = connect(**db_config)
    try:
        with connection.cursor() as cursor:
            cursor.execute(UPDATE_QUERY, params)
            if cursor.rowcount == 1:
                connection.commit()
                return "set", VALIDATOR_2
            # Nothing updated: either the row is absent or someone already
            # signed off. Say which, so the upload report can tell them apart.
            cursor.execute(EXISTING_QUERY, params)
            existing = cursor.fetchone()
        connection.commit()
        if existing is None:
            return "no_row", None
        return "preserved", existing[1]
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file")
    parser.add_argument("validation_record")
    args = parser.parse_args()
    try:
        record = read_record(args.validation_record)
        full_seqid = record["full_seqid"]
        if not record["submission_ready"]:
            print(f"ℹ️ {full_seqid} not submission_ready — skipping validator_2 write.")
            return 0
        status, detail = set_qc_validator(record, load_db_config(args.config_file))
        if status == "set":
            print(
                f"✅ Success: lca_validation validator_2 set to '{VALIDATOR_2}' "
                f"for {full_seqid}"
            )
        elif status == "preserved":
            print(
                f"⚠️ Existing validator_2 preserved for {full_seqid}: "
                f"validator_2='{detail}'"
            )
        else:
            print(f"⚠️ No lca_validation row for {full_seqid} — validator_2 not set.")
        return 0
    except Exception as error:
        print(f"❌ Database error: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
