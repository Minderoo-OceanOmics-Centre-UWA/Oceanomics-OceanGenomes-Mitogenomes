#!/usr/bin/env python3
"""Upsert the ``lca_validation`` row from the record written by species_validation.py.

This is the DB half of what used to be a single species_validation.py. Splitting it
out is what lets the species comparison -- which produces Found_in_blast_YN, and so
decides the whole QC gate -- run when uploads are skipped or no database is
configured at all. The comparison writes a JSON record; this reads it and writes the
row. Nothing here recomputes a verdict.

The record carries an ``action``:
  upsert  write the row described by the record
  skip    the comparison decided no row should be written (sample not validated, or
          no assembly key could be derived); this exits 0 having done nothing

Usage:
    push_species_validation.py [--force] <config.cfg> <validation_record.json>
"""

import argparse
import json
import sys

try:
    import psycopg2
except ImportError:  # Allows unit tests to inject a connection factory.
    psycopg2 = None

# load_db_config is identical across every push script in bin/; reuse rather
# than copy. bin/ is on PATH at runtime and sys.path[0] is this file's
# directory, so the sibling import resolves the same way species_validation.py
# imports species_name_utils.
from push_ena_validation_results import load_db_config


def upsert_lca_validation(
    db_params, key, validated_species_name, validator="nf-core", force=False,
    validated_rank=None, lca_genus=None
):
    """
    Insert/update the lca_validation row for this sample with the validated
    species name and validator tag. The composite key matches the assembly
    naming convention: (og_id, tech, seq_date, code, annotation).

    Guard: if an existing row carries a validator tag other than "nf-core",
    the row was probably set by a human reviewer or another pipeline and
    must not be silently overwritten by an automated re-run. Skip the
    upsert in that case unless ``force=True``.
    """
    og_id, tech, seq_date, code, annotation = key
    params = {
        "og_id": og_id,
        "tech": tech,
        "seq_date": seq_date,
        "code": code,
        "annotation": annotation,
        "validated_species_name": validated_species_name,
        "validator": validator,
        # The rank the evidence actually supported. A relaxed gate stays honest
        # only if what was relaxed is recorded, and a genus-level release must be
        # distinguishable from a species-level one after the fact.
        "validated_rank": validated_rank,
        # For a family-level release, the genus the LCA resolved. Submitting
        # 'Ophidiidae sp.' when the pipeline already resolved Lamprogrammus throws
        # away information, so record it: these surface as label-upgrade candidates
        # for a curator, and validated_rank makes them selectable if a later policy
        # wants to hold family matches for curation instead.
        "lca_genus": lca_genus,
    }
    conn = None
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            if not force:
                cur.execute(
                    """
                    SELECT validator, validated_species_name
                    FROM lca_validation
                    WHERE og_id = %(og_id)s
                      AND tech = %(tech)s
                      AND seq_date = %(seq_date)s
                      AND code = %(code)s
                      AND annotation = %(annotation)s
                    """,
                    params,
                )
                existing = cur.fetchone()
                if existing is not None:
                    existing_validator, existing_species = existing
                    if (
                        existing_validator is not None
                        and str(existing_validator).strip().lower() != "nf-core"
                    ):
                        print(
                            "⚠️ Existing values preserved for lca_validation "
                            f"{og_id}.{tech}.{seq_date}.{code}.{annotation}: "
                            f"validator='{existing_validator}' (not 'nf-core'); "
                            "pass --force to overwrite."
                        )
                        print(
                            "📌 Preserved stored values: "
                            f"{{'validator': '{existing_validator}', "
                            f"'validated_species_name': '{existing_species}'}}"
                        )
                        return True

            upsert_query = """
            INSERT INTO lca_validation (
                og_id, tech, seq_date, code, annotation,
                validated_species_name, validator, validated_rank, lca_genus
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s,
                %(validated_species_name)s, %(validator)s, %(validated_rank)s,
                %(lca_genus)s
            )
            ON CONFLICT (og_id, tech, seq_date, code, annotation)
            DO UPDATE SET
                validated_species_name = EXCLUDED.validated_species_name,
                validator = EXCLUDED.validator,
                validated_rank = EXCLUDED.validated_rank,
                lca_genus = EXCLUDED.lca_genus
            """
            cur.execute(upsert_query, params)
        conn.commit()
        if force:
            print(
                "✅ Success: lca_validation overwritten (--force) for "
                f"{og_id}.{tech}.{seq_date}.{code}.{annotation} "
                f"-> validated_species_name='{validated_species_name}', validator='{validator}'"
            )
        else:
            print(
                "✅ Success: lca_validation upserted for "
                f"{og_id}.{tech}.{seq_date}.{code}.{annotation} "
                f"-> validated_species_name='{validated_species_name}', validator='{validator}'"
            )
        return True
    except Exception as e:
        if conn is not None:
            conn.rollback()
        print(f"❌ Database error while writing lca_validation: {e}")
        return False
    finally:
        if conn is not None:
            conn.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite an existing lca_validation row even when its validator "
            "is something other than 'nf-core' (e.g. a human reviewer). Off "
            "by default."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("record_file")
    args = parser.parse_args()

    with open(args.record_file) as handle:
        record = json.load(handle)

    action = record.get("action")
    if action == "skip":
        print(
            f"ℹ️ No lca_validation row to write: {record.get('reason') or 'skipped'}"
        )
        return 0
    if action != "upsert":
        print(f"❌ Unknown action '{action}' in {args.record_file}")
        return 1

    key_fields = record.get("key")
    if not key_fields:
        print(f"❌ Record asks for an upsert but carries no key: {args.record_file}")
        return 1
    key = (
        key_fields["og_id"], key_fields["tech"], key_fields["seq_date"],
        key_fields["code"], key_fields["annotation"],
    )

    db_params = load_db_config(args.config_file)
    ok = upsert_lca_validation(
        db_params,
        key,
        record.get("validated_species_name"),
        validator=record.get("validator") or "nf-core",
        force=args.force,
        validated_rank=record.get("validated_rank"),
        lca_genus=record.get("lca_genus"),
    )
    # Non-zero on a failed DB write so Nextflow sees the failure. A failed
    # lca_validation write used to print and return normally, so the task exited 0
    # and Nextflow saw success -- which is exactly how rows lost to the foreign-key
    # write-ordering race stayed invisible until someone read a per-sample upload
    # log. The sibling pushers already fail loudly; this keeps the third writer in
    # line now that it lives in its own process.
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
