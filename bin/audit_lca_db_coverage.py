#!/usr/bin/env python3
"""Report published LCA results that never reached the database.

The upload scripts used to catch a per-row error, print a tally and exit 0 even
when PostgreSQL had aborted the transaction and rolled the whole batch back, so
a lost upload leaves no trace in the Nextflow trace: the task shows COMPLETED
and its published ``.upload.txt`` ends in a tick. The only reliable way to find
those losses is to compare what a run published on disk against what is in the
database.

For every assembly directory under ``<run>/mitogenomes/<OG>/<assembly>/`` this
compares the row counts in ``lca/`` against ``blast_filtered_lca``, ``lca`` and
``lca_raw_results``. An assembly is only reported when the published file has
data and the table has none: a sample whose regions produced no hits above
threshold writes empty ``blast.*.filtered.tsv`` and header-only ``lca_raw.*.tsv``
files, and correctly has no rows.

Usage:
    audit_lca_db_coverage.py <config.cfg> <run_dir> [<run_dir> ...]
    audit_lca_db_coverage.py --json gaps.json <config.cfg> /scratch/.../batch-*
"""

import argparse
import configparser
import json
import sys
from collections import Counter
from pathlib import Path

import psycopg2

# Table -> the published files that feed it, and whether those files have a
# header line that must not be counted as data.
SOURCES = {
    "blast_filtered_lca": ("blast.*.filtered.tsv", False),
    "lca": ("lca.*.tsv", True),
    "lca_raw_results": ("lca_raw.*.tsv", True),
}


def load_db_config(config_file):
    config = configparser.ConfigParser()
    if not config.read(config_file) or not config.has_section("postgres"):
        raise ValueError(f"Missing [postgres] configuration in {config_file}")
    return {
        "dbname": config.get("postgres", "dbname"),
        "user": config.get("postgres", "user"),
        "password": config.get("postgres", "password"),
        "host": config.get("postgres", "host"),
        "port": config.getint("postgres", "port"),
    }


def load_db_counts(db_params):
    """Row counts per assembly key, for each table this audit checks."""
    counts = {}
    with psycopg2.connect(**db_params) as conn:
        with conn.cursor() as cur:
            for table in SOURCES:
                cur.execute(
                    f"SELECT og_id, tech, seq_date, code, count(*) "
                    f"FROM {table} GROUP BY 1, 2, 3, 4"
                )
                counts[table] = {
                    (r[0], r[1], str(r[2]), r[3]): r[4] for r in cur.fetchall()
                }
    return counts


def count_data_rows(paths, has_header):
    """Data rows across a set of TSVs, ignoring headers and blank lines."""
    total = 0
    for path in paths:
        with open(path) as handle:
            lines = [line for line in handle if line.strip()]
        total += max(0, len(lines) - (1 if has_header else 0))
    return total


def assembly_key(assembly_dir):
    """('OG2343', 'ilmn', '260514', 'getorg1770') from the directory name."""
    parts = assembly_dir.name.split(".")
    return tuple(parts) if len(parts) == 4 else None


def find_gaps(run_dirs, db_counts):
    """Assemblies whose published LCA output is missing from the database."""
    gaps = []
    for run_dir in run_dirs:
        for assembly_dir in sorted(Path(run_dir).glob("mitogenomes/*/*/")):
            key = assembly_key(assembly_dir)
            if key is None or not (assembly_dir / "lca").is_dir():
                continue
            missing = {}
            for table, (pattern, has_header) in SOURCES.items():
                published = count_data_rows(
                    sorted((assembly_dir / "lca").glob(pattern)), has_header
                )
                if published and not db_counts[table].get(key, 0):
                    missing[table] = published
            if missing:
                gaps.append(
                    {
                        "run": Path(run_dir).name,
                        "assembly": assembly_dir.name,
                        "assembly_dir": str(assembly_dir),
                        "og_id": key[0],
                        "missing": missing,
                    }
                )
    return gaps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file")
    parser.add_argument("run_dirs", nargs="+", help="Pipeline outdirs to audit.")
    parser.add_argument("--json", help="Also write the gaps to this file.")
    args = parser.parse_args()

    gaps = find_gaps(args.run_dirs, load_db_counts(load_db_config(args.config_file)))

    for run in sorted({gap["run"] for gap in gaps}):
        rows = [gap for gap in gaps if gap["run"] == run]
        print(f"== {run} ({len(rows)})")
        for gap in rows:
            detail = "; ".join(
                f"{table} (published={n})" for table, n in sorted(gap["missing"].items())
            )
            print(f"   {gap['assembly']}\t{detail}")
    by_table = Counter(t for gap in gaps for t in gap["missing"])
    print(f"\nTOTAL {len(gaps)} assemblies with gaps", dict(by_table) if gaps else "")

    if args.json:
        Path(args.json).write_text(json.dumps(gaps, indent=1))
        print(f"wrote {args.json}")
    return 1 if gaps else 0


if __name__ == "__main__":
    sys.exit(main())
