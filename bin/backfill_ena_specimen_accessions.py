#!/usr/bin/env python3
"""Populate the ENA specimen registry from existing OceanOmics specimen data."""

from __future__ import annotations

import argparse
import configparser
import re
import sys
from collections import defaultdict
from pathlib import Path

import psycopg2


def load_config(path: Path) -> dict[str, object]:
    parser = configparser.ConfigParser()
    parser.read(path)
    return {
        "dbname": parser.get("postgres", "dbname"),
        "user": parser.get("postgres", "user"),
        "password": parser.get("postgres", "password"),
        "host": parser.get("postgres", "host"),
        "port": parser.getint("postgres", "port"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    try:
        with psycopg2.connect(**load_config(Path(args.config))) as connection:
            with connection.cursor() as cursor:
                cursor.execute(
                    """
                    SELECT DISTINCT og_id FROM sample
                    UNION
                    SELECT DISTINCT og_id FROM mitogenome_data
                    """
                )
                og_ids = sorted(
                    row[0]
                    for row in cursor.fetchall()
                    if row[0] and re.fullmatch(r"OG[0-9]+", row[0])
                )
                numeric_to_og: dict[int, str] = {}
                for og_id in og_ids:
                    numeric = int(og_id[2:])
                    if numeric > 999999:
                        # Matches the ena_specimen_og_numeric_check constraint in
                        # sql/004_ena_candidate_packages.sql, so an over-wide OG
                        # fails here with a name rather than as a check violation.
                        raise ValueError(
                            f"{og_id} exceeds the six-digit og_numeric range"
                        )
                    if numeric in numeric_to_og and numeric_to_og[numeric] != og_id:
                        raise ValueError(
                            f"OG numeric collision: {numeric_to_og[numeric]} and {og_id}"
                        )
                    numeric_to_og[numeric] = og_id
                # BioSample is one INSDC namespace; the prefix records only which
                # archive minted the accession (SAMEA/EBI, SAMN/NCBI, SAMD/DDBJ).
                # OceanOmics registers at NCBI, so a SAMEA-only filter matches
                # nothing.
                cursor.execute(
                    """
                    SELECT og_id, ncbi_biosample_id
                    FROM sample
                    WHERE ncbi_biosample_id ~ '^SAM(EA|N|D)[0-9]+$'
                    """
                )
                sample_accessions = {og_id: accession for og_id, accession in cursor.fetchall()}
                cursor.execute(
                    """
                    SELECT og_id, biosample_accession
                    FROM draft_genomes
                    WHERE biosample_accession ~ '^SAM(EA|N|D)[0-9]+$'
                    """
                )
                draft_accessions: dict[str, set[str]] = defaultdict(set)
                for og_id, accession in cursor.fetchall():
                    draft_accessions[og_id].add(accession)
                # sample is the source of truth, so a draft_genomes disagreement
                # is reported and overridden rather than aborting the backfill.
                accession_sets: dict[str, set[str]] = defaultdict(set)
                for og_id, accession in sample_accessions.items():
                    accession_sets[og_id].add(accession)
                for og_id, values in draft_accessions.items():
                    authoritative = sample_accessions.get(og_id)
                    if authoritative is None:
                        if len(values) > 1:
                            print(
                                f"WARNING: conflicting draft_genomes biosample "
                                f"og_id={og_id} values={sorted(values)}; skipping",
                                file=sys.stderr,
                            )
                            continue
                        accession_sets[og_id].update(values)
                    elif values != {authoritative}:
                        print(
                            f"WARNING: conflicting biosample og_id={og_id} "
                            f"sample={authoritative} draft_genomes={sorted(values)}; "
                            f"using sample",
                            file=sys.stderr,
                        )
                print(
                    f"specimens={len(og_ids)} "
                    f"ena_biosamples={len(accession_sets)} apply={args.apply}"
                )
                for og_id, values in sorted(accession_sets.items()):
                    print(f"biosample={og_id}:{next(iter(values))}")
                if args.apply:
                    for og_id in og_ids:
                        accession = next(iter(accession_sets.get(og_id, [])), None)
                        cursor.execute(
                            """
                            INSERT INTO ena_specimen_accessions (
                                og_id, og_numeric, ena_biosample_accession,
                                accession_source, verified_at
                            ) VALUES (%s, %s, %s, %s, CURRENT_TIMESTAMP)
                            ON CONFLICT (og_id) DO UPDATE SET
                                ena_biosample_accession = COALESCE(
                                    EXCLUDED.ena_biosample_accession,
                                    ena_specimen_accessions.ena_biosample_accession
                                ),
                                accession_source = COALESCE(
                                    EXCLUDED.accession_source,
                                    ena_specimen_accessions.accession_source
                                ),
                                verified_at = CASE
                                    WHEN EXCLUDED.ena_biosample_accession IS NOT NULL
                                    THEN CURRENT_TIMESTAMP
                                    ELSE ena_specimen_accessions.verified_at
                                END,
                                updated_at = CURRENT_TIMESTAMP
                            """,
                            (
                                og_id,
                                int(og_id[2:]),
                                accession,
                                (
                                    "sample"
                                    if og_id in sample_accessions
                                    else "draft_genomes"
                                )
                                if accession
                                else None,
                            ),
                        )
            if not args.apply:
                connection.rollback()
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
