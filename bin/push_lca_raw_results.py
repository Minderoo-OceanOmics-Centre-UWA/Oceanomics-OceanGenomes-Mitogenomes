#!/usr/bin/env python3
"""
Push raw LCA / BLAST hit rows into the OceanOmics ``lca_raw_results`` table.

The LCA module emits one ``lca_raw.<region>.<og_id>.<tech>.<seq_date>.<code>.<annotation>.tsv``
file per (sample, region). Each row in the TSV represents one BLAST hit that
contributed to the LCA decision for that region.

Rows are keyed by their content: the ``lca_raw_results_content_unique``
constraint covers ``og_id, tech, seq_date, code, annotation, sequence_region,
accession_id`` plus a ``content_hash`` maintained by a DB trigger. So a re-run
that reproduces a hit exactly just refreshes that row's ``lca_run_date``, while a
re-run that produces a different result for the same hit is inserted as a new row
and the earlier one is kept.

``--force`` additionally discards that history: after uploading, rows for the
same sample/region that this run did not produce are deleted, leaving only the
current result. The ``og_num`` and ``content_hash`` columns are maintained by the
DB and are intentionally not written here.

Usage:
    push_lca_raw_results.py [--force] <config.cfg> <sample> <lca_raw.tsv> [<lca_raw.tsv> ...]
"""

import argparse
import configparser
import sys
from pathlib import Path

import psycopg2
import pandas as pd
import numpy as np

from pg_row_guard import row_savepoint


# Map TSV header -> DB column. Anything not listed here is dropped.
# Note: the DB column is `specific_epiphet` (typo retained from the schema);
# the TSV header is the camelCase `specificEpithet`.
TSV_TO_DB = {
    "domain":                       "domain",
    "phylum":                       "phylum",
    "class":                        "class",
    "order":                        "order",
    "family":                       "family",
    "genus":                        "genus",
    "specificEpithet":              "specific_epiphet",
    "scientificName":               "scientific_name",
    "scientificNameAuthorship":     "scientific_name_authorship",
    "taxonRank":                    "taxon_rank",
    "taxonID":                      "taxon_id",
    "taxonID_db":                   "taxon_id_db",
    "verbatimIdentification":       "verbatim_identification",
    "accession_id":                 "accession_id",
    "accession_id_ref_db":          "accession_id_ref_db",
    "taxonRank_db":                 "taxon_rank_db",
    "percent_match":                "percent_match",
    "percent_query_cover":          "percent_query_cover",
    "percent_query_cover_hsp":      "percent_query_cover_hsp",
    "alignment_length":             "alignment_length",
    "subject_length":               "subject_length",
    "sequence_length":              "sequence_length",
    "confidence_score":             "confidence_score",
    "sequence_region":              "sequence_region",
    "lca_run_date":                 "lca_run_date",
}

# Columns supplied per row in the INSERT, in declaration order.
KEY_COLS = ["og_id", "tech", "seq_date", "code", "annotation"]
NON_KEY_DB_COLS = [
    "sequence_region", "lca_run_date",
    "domain", "phylum", "class", "order", "family", "genus",
    "specific_epiphet", "scientific_name", "scientific_name_authorship",
    "taxon_rank", "taxon_id", "taxon_id_db", "verbatim_identification",
    "accession_id", "accession_id_ref_db", "taxon_rank_db",
    "percent_match", "percent_query_cover", "percent_query_cover_hsp",
    "alignment_length", "subject_length", "sequence_length",
    "confidence_score",
]
ALL_DB_COLS = KEY_COLS + NON_KEY_DB_COLS

# Conflict target for the upsert. Named rather than inferred because part of
# the constraint (content_hash) is populated by a BEFORE INSERT trigger.
CONFLICT_CONSTRAINT = "lca_raw_results_content_unique"

# The columns that constraint covers, alongside the trigger-maintained
# content_hash. lca_run_date is deliberately NOT part of it: a re-run with an
# identical result must update the existing row's date rather than insert a
# second copy of it.
UNIQUE_KEY = [
    "og_id", "tech", "seq_date", "code", "annotation",
    "sequence_region", "accession_id",
]

# The group --force prunes: everything this sample/region has ever recorded.
HISTORY_KEY = [
    "og_id", "tech", "seq_date", "code", "annotation", "sequence_region",
]


def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    cfg = configparser.ConfigParser()
    cfg.read(config_file)
    if not cfg.has_section("postgres"):
        raise ValueError("❌ Missing [postgres] section in config file.")
    return {
        "dbname":   cfg.get("postgres", "dbname"),
        "user":     cfg.get("postgres", "user"),
        "password": cfg.get("postgres", "password"),
        "host":     cfg.get("postgres", "host"),
        "port":     cfg.getint("postgres", "port"),
    }


def parse_assembly_key(seq_id):
    """og_id.tech.seq_date.code.annotation -> (og_id, tech, seq_date, code, annotation)."""
    parts = str(seq_id).split(".")
    if len(parts) < 5:
        return None
    return tuple(parts[:5])


def normalise(v):
    """Coerce pandas NaN/empty to None so psycopg2 binds them as NULL."""
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    if isinstance(v, str) and v.strip() == "":
        return None
    return v


def build_row_params(tsv_row):
    """Map one TSV row dict to the parameter dict required by the INSERT."""
    seq_id = tsv_row.get("seq_id")
    key = parse_assembly_key(seq_id)
    if key is None:
        return None
    og_id, tech, seq_date, code, annotation = key

    params = {
        "og_id":     og_id,
        "tech":      tech,
        "seq_date":  seq_date,
        "code":      code,
        "annotation": annotation,
    }
    for tsv_col, db_col in TSV_TO_DB.items():
        params[db_col] = normalise(tsv_row.get(tsv_col))
    return params


def build_insert_query():
    # Quote every column name. Several columns (`order`, `class`) collide
    # with reserved SQL keywords, and quoting them all keeps the script
    # robust to any future column additions that might do the same.
    def q(col):
        return f'"{col}"'

    cols_sql = ", ".join(q(c) for c in ALL_DB_COLS)
    placeholders = ", ".join(f"%({c})s" for c in ALL_DB_COLS)

    # A conflict means this run reproduced an existing hit exactly, so the only
    # new information is that it was seen again. A hit whose values changed
    # hashes differently, misses the constraint, and lands as a new row.
    # RETURNING (xmax = 0) distinguishes the two outcomes: xmax is 0 on a fresh
    # insert and non-zero on the DO UPDATE path.
    return f"""
    INSERT INTO lca_raw_results ({cols_sql})
    VALUES ({placeholders})
    ON CONFLICT ON CONSTRAINT {CONFLICT_CONSTRAINT}
    DO UPDATE SET "lca_run_date" = EXCLUDED."lca_run_date"
    RETURNING (xmax = 0) AS inserted
    """


def prune_superseded(cur, written):
    """Delete rows this run did not produce, for each sample/region it touched.

    `written` maps a HISTORY_KEY tuple to the set of lca_run_date values this run
    wrote for it. Identical-content rows have just had their date refreshed to the
    current run's, so anything left on another date is a superseded result.
    """
    where = " AND ".join(f'"{c}" = %s' for c in HISTORY_KEY)
    query = (
        f"DELETE FROM lca_raw_results WHERE {where} "
        f"AND \"lca_run_date\"::text <> ALL(%s)"
    )
    deleted = 0
    for key, run_dates in written.items():
        cur.execute(query, list(key) + [[str(d) for d in run_dates]])
        deleted += cur.rowcount
    return deleted


def process_file(cur, path, written):
    """Upload one lca_raw TSV. Records what it wrote into `written` for --force."""
    print(f"📂 Reading lca_raw file: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str).replace({np.nan: None})
    if df.empty:
        print(f"ℹ️ {path} has no rows; nothing to insert.")
        return 0, 0, 0

    insert_query = build_insert_query()
    inserted = refreshed = failed = 0

    for _, row in df.iterrows():
        params = build_row_params(row.to_dict())
        if params is None:
            failed += 1
            print(f"⚠️ Skipped row with unparseable seq_id: {row.get('seq_id')!r}")
            continue
        try:
            # The RETURNING row has to be read before the savepoint is released:
            # RELEASE SAVEPOINT is itself a statement on this cursor and would
            # replace the INSERT's result set.
            with row_savepoint(cur):
                cur.execute(insert_query, params)
                was_insert = cur.fetchone()[0]
            if was_insert:
                inserted += 1
            else:
                refreshed += 1
            key = tuple(params[c] for c in HISTORY_KEY)
            written.setdefault(key, set()).add(params["lca_run_date"])
        except Exception as e:
            failed += 1
            print(
                f"❌ Insert failed for accession {params.get('accession_id')} "
                f"in region {params.get('sequence_region')}: {e}"
            )

    return inserted, refreshed, failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "After uploading, delete rows for the same sample/region that this "
            "run did not produce, so only the current result is kept. Off by "
            "default: superseded results are preserved alongside the new ones."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("sample", help="OG id of the sample being processed (for logging).")
    parser.add_argument("lca_raw_files", nargs="+", help="One or more lca_raw.*.tsv files.")
    args = parser.parse_args()

    db_params = load_db_config(args.config_file)

    total_inserted = total_refreshed = total_failed = 0
    total_pruned = 0
    # HISTORY_KEY tuple -> set of lca_run_date values written for it, across all
    # input files, so a region split over several files prunes correctly.
    written = {}
    conn = None
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            for f in args.lca_raw_files:
                ins, refr, fld = process_file(cur, f, written)
                total_inserted += ins
                total_refreshed += refr
                total_failed += fld
            # Only prune once every file has been uploaded, and only if the
            # upload was clean: pruning around a failed row could delete the
            # older result without a replacement having landed.
            if args.force and total_failed == 0:
                total_pruned = prune_superseded(cur, written)
            elif args.force:
                print("⚠️ Skipping --force prune: upload had failures.")
        conn.commit()
    except Exception as e:
        if conn is not None:
            conn.rollback()
        print(f"❌ Database error: {e}")
        total_failed += 1
    finally:
        if conn is not None:
            conn.close()

    mode = "replace-history" if args.force else "keep-history"
    tally = f"{total_inserted} inserted, {total_refreshed} refreshed"
    if args.force:
        tally += f", {total_pruned} superseded rows deleted"

    if total_failed == 0:
        print(
            f"✅ Success: lca_raw_results upload complete for {args.sample} "
            f"({mode}): {tally}."
        )
        return 0

    print(
        f"⚠️ lca_raw_results upload finished with errors for {args.sample} "
        f"({mode}): {tally}, {total_failed} failed."
    )
    # Non-zero so Nextflow surfaces a partial upload instead of publishing a
    # green task whose log quietly reports missing rows.
    return 1


if __name__ == "__main__":
    sys.exit(main())
