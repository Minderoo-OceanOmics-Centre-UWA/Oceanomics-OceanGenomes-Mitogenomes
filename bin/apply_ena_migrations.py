#!/usr/bin/env python3
"""Audit and apply the ordered OceanOmics ENA PostgreSQL migrations."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import sys
from pathlib import Path

import psycopg2


MIGRATIONS = [
    "001_create_ena_validation_attempts.sql",
    "002_ena_validation_attempts_single_row_per_attempt.sql",
    "003_mitogenome_data_uniform_depth.sql",
    "004_ena_candidate_packages.sql",
    "005_ena_validation_attempts_genome_context.sql",
    "006_ena_tech_aware_locus_tags.sql",
    "007_insdc_biosample_accessions.sql",
    "008_ena_submission_queue.sql",
    "009_ena_locus_registry_canonical_order.sql",
    "010_drop_ena_locus_tables.sql",
    "011_drop_local_package_validation.sql",
    "012_drop_ena_selection_layer.sql",
]


def load_config(path: Path) -> dict[str, object]:
    parser = configparser.ConfigParser()
    if not parser.read(path) or not parser.has_section("postgres"):
        raise ValueError(f"Missing [postgres] configuration in {path}")
    return {
        "dbname": parser.get("postgres", "dbname"),
        "user": parser.get("postgres", "user"),
        "password": parser.get("postgres", "password"),
        "host": parser.get("postgres", "host"),
        "port": parser.getint("postgres", "port"),
    }


def audit(cursor) -> dict[str, object]:
    cursor.execute(
        """
        SELECT
            to_regclass('public.ena_validation_attempts') IS NOT NULL,
            to_regclass('public.ena_candidate_packages') IS NULL
                AND to_regclass('public.ena_submission_selections') IS NULL
                AND to_regclass('public.ena_submission_queue') IS NULL,
            to_regclass('public.ena_locus_registry') IS NULL
                AND to_regclass('public.ena_candidate_loci') IS NULL
        """
    )
    tables = cursor.fetchone()
    cursor.execute(
        """
        SELECT EXISTS (
            SELECT 1 FROM information_schema.columns
            WHERE table_schema = 'public'
              AND table_name = 'mitogenome_data'
              AND column_name = 'mean_depth'
        )
        """
    )
    mean_depth = cursor.fetchone()[0]
    measured_depth_rows = 0
    if mean_depth:
        cursor.execute(
            """
            SELECT count(*)
            FROM mitogenome_data
            WHERE mean_depth IS NOT NULL
            """
        )
        measured_depth_rows = cursor.fetchone()[0]
    return {
        "validation_table": tables[0],
        # Migration 012 retires the selection layer: choosing and submitting a
        # package belongs to the downstream submission pipeline, so the absence
        # of these three is the healthy state.
        "selection_tables_dropped": tables[1],
        # Migration 010 retires the locus tables: tag allocation belongs to the
        # downstream submission pipeline, so their absence is the healthy state.
        "locus_tables_dropped": tables[2],
        "mean_depth_column": mean_depth,
        "measured_depth_rows": measured_depth_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument(
        "--migration-dir",
        default=str(Path(__file__).resolve().parents[1] / "sql"),
    )
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    try:
        migration_dir = Path(args.migration_dir)
        migrations = [migration_dir / name for name in MIGRATIONS]
        missing = [str(path) for path in migrations if not path.is_file()]
        if missing:
            raise ValueError(f"Missing migration files: {', '.join(missing)}")
        with psycopg2.connect(**load_config(Path(args.config))) as connection:
            connection.autocommit = True
            with connection.cursor() as cursor:
                cursor.execute("SELECT pg_advisory_lock(hashtext(%s))", ("oceanomics_ena_migrations",))
                try:
                    try:
                        before = audit(cursor)
                        print(f"before={before}")
                        if not args.check_only:
                            for path in migrations:
                                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                                print(f"applying={path.name} sha256={digest}")
                                cursor.execute(path.read_text())
                            after = audit(cursor)
                            print(f"after={after}")
                            if not all(
                                after[key]
                                for key in (
                                    "validation_table",
                                    "selection_tables_dropped",
                                    "locus_tables_dropped",
                                    "mean_depth_column",
                                )
                            ):
                                raise RuntimeError("Post-migration ENA schema audit failed")
                        else:
                            print("check_only=true")
                    except Exception:
                        cursor.execute("ROLLBACK")
                        raise
                finally:
                    cursor.execute(
                        "SELECT pg_advisory_unlock(hashtext(%s))",
                        ("oceanomics_ena_migrations",),
                    )
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
