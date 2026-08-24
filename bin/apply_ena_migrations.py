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
    "013_drop_ena_candidate_runs.sql",
    "014_ena_validation_attempts_full_seqid.sql",
    "015_ena_submissions.sql",
    "016_mitogenome_data_og_num_generated.sql",
    "017_lca_content_addressed_rows.sql",
    "018_mitogenome_data_og_num_first.sql",
    "019_ena_validation_attempts_recompute_submission_ready.sql",
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
                AND to_regclass('public.ena_submission_queue') IS NULL
                AND to_regclass('public.ena_candidate_runs') IS NULL,
            to_regclass('public.ena_locus_registry') IS NULL
                AND to_regclass('public.ena_candidate_loci') IS NULL,
            to_regclass('public.ena_submissions') IS NOT NULL
                AND to_regclass('public.ena_submission_status') IS NOT NULL
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
    cursor.execute(
        """
        SELECT EXISTS (
            SELECT 1 FROM pg_attribute
            WHERE attrelid = to_regclass('public.mitogenome_data')
              AND attname = 'og_num'
              AND attgenerated = 's'
              AND attnum = 1
        )
        """
    )
    og_num_generated = cursor.fetchone()[0]
    cursor.execute(
        """
        SELECT
            (SELECT count(*) FROM information_schema.columns
              WHERE table_schema = 'public'
                AND table_name IN ('lca', 'lca_raw_results')
                AND column_name = 'content_hash'
                AND is_nullable = 'NO') = 2
            AND EXISTS (
                SELECT 1 FROM pg_proc p JOIN pg_namespace n ON n.oid = p.pronamespace
                WHERE n.nspname = 'public' AND p.proname = 'lca_set_content_hash')
            AND (SELECT count(*) FROM pg_constraint
                  WHERE conname IN ('lca_content_unique',
                                    'lca_raw_results_content_unique')) = 2
        """
    )
    lca_content_addressing = cursor.fetchone()[0]
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
        # Migrations 012 and 013 retire the selection layer: choosing and
        # submitting a package, and the run accessions that identify the
        # raw-read submissions, belong to the downstream submission pipeline,
        # so the absence of all four is the healthy state.
        "selection_tables_dropped": tables[1],
        # Migration 010 retires the locus tables: tag allocation belongs to the
        # downstream submission pipeline, so their absence is the healthy state.
        "locus_tables_dropped": tables[2],
        # Migration 015 adds the submission ledger and the status view. This
        # pipeline never writes them; the downstream submitter does, and it can
        # only do that if they exist.
        "submission_ledger": tables[3],
        "mean_depth_column": mean_depth,
        "measured_depth_rows": measured_depth_rows,
        # Migration 016 restores og_num as a stored generated column, matching
        # every other table in the schema. Nothing writes it, so if the
        # generation expression is missing the column silently fills with NULLs.
        # Migration 018 additionally pins it to column 1, where it is meant to
        # be: a rebuild that moves it is the same accident that lost the
        # generation expression the first time, so the audit checks both.
        "og_num_generated": og_num_generated,
        # Migration 017 records the content-addressing rework of lca and
        # lca_raw_results. The push scripts name lca_content_unique and
        # lca_raw_results_content_unique as ON CONFLICT targets, so without
        # these the first LCA upload of a run fails outright.
        "lca_content_addressing": lca_content_addressing,
    }


def ensure_ledger(cursor) -> None:
    """Create the applied-migrations ledger if it is not there yet.

    This cannot live in sql/ as a migration of its own: it has to exist before
    any migration can be recorded as applied.
    """
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS public.schema_migrations (
            filename   TEXT PRIMARY KEY,
            sha256     TEXT NOT NULL,
            applied_at TIMESTAMPTZ NOT NULL DEFAULT now()
        )
        """
    )


def read_ledger(cursor) -> dict[str, str]:
    cursor.execute("SELECT filename, sha256 FROM public.schema_migrations")
    return dict(cursor.fetchall())


def check_drift(migrations, applied: dict[str, str], force: bool) -> None:
    """Refuse to run when an already-applied migration file has been edited.

    The database was built from the bytes recorded in the ledger, so a file that
    has changed since no longer describes the schema it produced. Editing an
    applied migration means writing a new one instead.
    """
    drifted = [
        (path.name, applied[path.name], hashlib.sha256(path.read_bytes()).hexdigest())
        for path in migrations
        if path.name in applied
        and applied[path.name] != hashlib.sha256(path.read_bytes()).hexdigest()
    ]
    for name, recorded, current in drifted:
        message = f"{name} was modified after it was applied: recorded={recorded} current={current}"
        if force:
            print(f"WARNING: {message}")
        else:
            raise RuntimeError(message + " (pass --force to apply anyway)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument(
        "--migration-dir",
        default=str(Path(__file__).resolve().parents[1] / "sql"),
    )
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument(
        "--baseline",
        action="store_true",
        help="Record every listed migration as applied WITHOUT running any of it. "
        "Used once, to onboard a database the migrations were already applied to.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Downgrade the modified-after-applied check to a warning.",
    )
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
                        ensure_ledger(cursor)
                        applied = read_ledger(cursor)

                        if args.baseline:
                            if applied:
                                raise RuntimeError(
                                    f"Ledger already records {len(applied)} migration(s); "
                                    "--baseline is only for onboarding an empty ledger"
                                )
                            for path in migrations:
                                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                                print(f"baselining={path.name} sha256={digest}")
                                cursor.execute(
                                    "INSERT INTO public.schema_migrations (filename, sha256) VALUES (%s, %s)",
                                    (path.name, digest),
                                )
                            print(f"baselined={len(migrations)} applied_nothing=true")
                            return 0

                        check_drift(migrations, applied, args.force)
                        pending = [path for path in migrations if path.name not in applied]

                        before = audit(cursor)
                        print(f"before={before}")
                        print(f"applied={len(applied)} pending={[path.name for path in pending]}")
                        if not args.check_only:
                            for path in migrations:
                                if path.name in applied:
                                    print(f"skipping={path.name} already_applied=true")
                                    continue
                                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                                print(f"applying={path.name} sha256={digest}")
                                cursor.execute(path.read_text())
                                cursor.execute(
                                    "INSERT INTO public.schema_migrations (filename, sha256) VALUES (%s, %s)",
                                    (path.name, digest),
                                )
                            after = audit(cursor)
                            print(f"after={after}")
                            if not all(
                                after[key]
                                for key in (
                                    "validation_table",
                                    "selection_tables_dropped",
                                    "locus_tables_dropped",
                                    "submission_ledger",
                                    "mean_depth_column",
                                    "og_num_generated",
                                    "lca_content_addressing",
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
