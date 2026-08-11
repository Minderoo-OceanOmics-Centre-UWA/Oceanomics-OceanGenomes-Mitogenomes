#!/usr/bin/env python3
"""Report non-secret database prerequisites for an ENA genome-package run."""

from __future__ import annotations

import argparse
import configparser
import csv
import re
import sys
from pathlib import Path

import psycopg2
from psycopg2 import sql


INSDC_BIOSAMPLE_PATTERN = re.compile(r"SAM(?:EA|N|D)[0-9]+")
ENA_BIOSAMPLE_PATTERN = re.compile(r"SAMEA[0-9]+")


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


def coverage_status(accession: str) -> str:
    """Bucket a specimen by what has to happen before ENA will accept it.

    Deliberately not a reachability check against the ENA browser: that API
    serves the EBI BioSamples mirror and answers 200 for NCBI-registered
    specimens that webin then refuses, so it would report every SAMN as ready.
    Only a SAMEA is resolvable by webin's submission sample service.
    """
    if not accession:
        return "MISSING"
    if ENA_BIOSAMPLE_PATTERN.fullmatch(accession):
        return "SUBMITTABLE"
    if INSDC_BIOSAMPLE_PATTERN.fullmatch(accession):
        return "REGISTERED_NCBI_ONLY"
    return "MALFORMED"


def write_coverage(cursor, path: Path) -> None:
    """Emit per-specimen BioSample readiness for every assembled mitogenome."""
    cursor.execute(
        """
        SELECT DISTINCT m.og_id, s.ncbi_biosample_id
        FROM mitogenome_data m
        LEFT JOIN sample s ON s.og_id = m.og_id
        ORDER BY m.og_id
        """
    )
    rows = cursor.fetchall()
    counts: dict[str, int] = {}
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["og_id", "biosample", "prefix", "status"])
        for og_id, raw in rows:
            accession = (raw or "").strip()
            status = coverage_status(accession)
            match = INSDC_BIOSAMPLE_PATTERN.match(accession) if accession else None
            prefix = re.match(r"SAM(?:EA|N|D)", accession).group(0) if match else ""
            counts[status] = counts.get(status, 0) + 1
            writer.writerow([og_id, accession, prefix, status])
    summary = " ".join(f"{key.lower()}={value}" for key, value in sorted(counts.items()))
    print(f"coverage_rows={len(rows)} {summary} coverage_tsv={path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--og-id")
    parser.add_argument(
        "--coverage-tsv",
        type=Path,
        help="write per-specimen BioSample readiness for every assembled mitogenome",
    )
    args = parser.parse_args()
    try:
        with psycopg2.connect(**load_config(Path(args.config))) as connection:
            with connection.cursor() as cursor:
                cursor.execute(
                    """
                    SELECT count(*), count(mean_depth),
                           count(*) FILTER (WHERE depth_method = 'remap_full_v1')
                    FROM mitogenome_data
                    """
                )
                total, measured, uniform = cursor.fetchone()
                print(
                    f"mitogenome_rows={total} mean_depth_rows={measured} "
                    f"uniform_depth_rows={uniform}"
                )
                cursor.execute(
                    """
                    SELECT count(*), count(ena_biosample_accession)
                    FROM ena_specimen_accessions
                    """
                )
                specimens, biosamples = cursor.fetchone()
                print(
                    f"ena_specimen_rows={specimens} "
                    f"ena_biosample_rows={biosamples}"
                )
                cursor.execute("SELECT count(*) FROM ena_locus_registry")
                print(f"ena_locus_rows={cursor.fetchone()[0]}")
                cursor.execute("SELECT count(*) FROM ena_candidate_packages")
                print(f"ena_candidate_rows={cursor.fetchone()[0]}")
                cursor.execute(
                    """
                    SELECT table_name, column_name
                    FROM information_schema.columns
                    WHERE table_schema = 'public'
                      AND column_name ILIKE '%biosample%'
                    ORDER BY table_name, column_name
                    """
                )
                sources = cursor.fetchall()
                print(
                    "biosample_columns="
                    + ",".join(f"{table}.{column}" for table, column in sources)
                )
                for table, column in sources:
                    query = sql.SQL(
                        """
                        SELECT count(*) FILTER (
                                   WHERE {column}::text ~ '^SAM(EA|N|D)[0-9]+$'
                               ),
                               count(*) FILTER (
                                   WHERE {column} IS NOT NULL
                                     AND btrim({column}::text) <> ''
                               )
                        FROM {table}
                        """
                    ).format(
                        column=sql.Identifier(column),
                        table=sql.Identifier(table),
                    )
                    cursor.execute(query)
                    insdc, populated = cursor.fetchone()
                    print(
                        f"biosample_source={table}.{column} "
                        f"insdc={insdc} populated={populated}"
                    )
                    if insdc:
                        cursor.execute(
                            """
                            SELECT EXISTS (
                                SELECT 1 FROM information_schema.columns
                                WHERE table_schema = 'public'
                                  AND table_name = %s
                                  AND column_name = 'og_id'
                            )
                            """,
                            (table,),
                        )
                        if cursor.fetchone()[0]:
                            query = sql.SQL(
                                """
                                SELECT DISTINCT og_id, {column}::text
                                FROM {table}
                                WHERE {column}::text ~ '^SAM(EA|N|D)[0-9]+$'
                                ORDER BY og_id
                                LIMIT 20
                                """
                            ).format(
                                column=sql.Identifier(column),
                                table=sql.Identifier(table),
                            )
                            cursor.execute(query)
                            print(
                                f"insdc_rows={table}."
                                + ",".join(
                                    f"{og_id}:{accession}"
                                    for og_id, accession in cursor.fetchall()
                                )
                            )
                if args.coverage_tsv:
                    write_coverage(cursor, args.coverage_tsv)
                if args.og_id:
                    cursor.execute(
                        """
                        SELECT s.nominal_species_id, sp.class, sp.ordr
                        FROM sample s
                        LEFT JOIN species sp
                          ON lower(sp.species) =
                             lower(trim(s.nominal_species_id))
                        WHERE s.og_id = %s
                        ORDER BY sp.ncbi_taxon_id NULLS LAST
                        LIMIT 1
                        """,
                        (args.og_id,),
                    )
                    taxonomy = cursor.fetchone()
                    if taxonomy:
                        print(
                            "taxonomy="
                            + "|".join(
                                "" if value is None else str(value)
                                for value in taxonomy
                            )
                        )
                    for table, column in sources:
                        cursor.execute(
                            """
                            SELECT EXISTS (
                                SELECT 1 FROM information_schema.columns
                                WHERE table_schema = 'public'
                                  AND table_name = %s
                                  AND column_name = 'og_id'
                            )
                            """,
                            (table,),
                        )
                        if cursor.fetchone()[0]:
                            query = sql.SQL(
                                """
                                SELECT DISTINCT {column}::text
                                FROM {table}
                                WHERE og_id = %s AND {column} IS NOT NULL
                                ORDER BY 1
                                """
                            ).format(
                                column=sql.Identifier(column),
                                table=sql.Identifier(table),
                            )
                            cursor.execute(query, (args.og_id,))
                            values = [row[0] for row in cursor.fetchall()]
                            if values:
                                print(
                                    f"og_biosample_source={table}.{column} "
                                    f"values={','.join(values)}"
                                )
                    cursor.execute(
                        """
                        SELECT tech, seq_date, code, mean_depth, depth_method
                        FROM mitogenome_data
                        WHERE og_id = %s
                        ORDER BY seq_date DESC, tech, code
                        """,
                        (args.og_id,),
                    )
                    for row in cursor.fetchall():
                        print(
                            "candidate="
                            + "|".join("" if value is None else str(value) for value in row)
                        )
                    cursor.execute(
                        """
                        SELECT table_name, column_name
                        FROM information_schema.columns
                        WHERE table_schema = 'public'
                          AND (
                              column_name ILIKE '%fastq%'
                              OR column_name ILIKE '%file_path%'
                              OR column_name ILIKE '%read_path%'
                          )
                          AND EXISTS (
                              SELECT 1
                              FROM information_schema.columns keyed
                              WHERE keyed.table_schema = 'public'
                                AND keyed.table_name =
                                    information_schema.columns.table_name
                                AND keyed.column_name = 'og_id'
                          )
                        ORDER BY table_name, column_name
                        """
                    )
                    for table, column in cursor.fetchall():
                        query = sql.SQL(
                            """
                            SELECT DISTINCT {column}::text
                            FROM {table}
                            WHERE og_id = %s AND {column} IS NOT NULL
                              AND btrim({column}::text) <> ''
                            ORDER BY 1
                            LIMIT 20
                            """
                        ).format(
                            column=sql.Identifier(column),
                            table=sql.Identifier(table),
                        )
                        cursor.execute(query, (args.og_id,))
                        values = [row[0] for row in cursor.fetchall()]
                        if values:
                            print(
                                f"og_input_source={table}.{column} "
                                f"values={','.join(values)}"
                            )
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
