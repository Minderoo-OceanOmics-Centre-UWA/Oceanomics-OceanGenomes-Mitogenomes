#!/usr/bin/env python3

import argparse
import psycopg2
import pandas as pd
import numpy as np
import configparser
import sys
import re
from pathlib import Path

# -------------------------------
# Load DB credentials from .cfg
# -------------------------------
def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    
    config = configparser.ConfigParser()
    config.read(config_file)

    if not config.has_section('postgres'):
        raise ValueError("❌ Missing [postgres] section in config file.")

    required_keys = ['dbname', 'user', 'password', 'host', 'port']
    for key in required_keys:
        if not config.has_option('postgres', key):
            raise ValueError(f"❌ Missing '{key}' in [postgres] section of config file.")

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')    
    }

# -------------------------------
# Helpers
# -------------------------------
def coerce_bool(val):
    """Coerce a TSV cell into True/False/None safely."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    return None

def try_parse_contig_stats(tsv_path: Path, target_contig="final_mitogenome"):
    """
    If tsv_path looks like contig_stats (has contig_id & was_circular),
    return stats and optional average coverage from the target row.
    If file missing or not a contig_stats file, return (None, None).
    """
    try:
        if not tsv_path.exists():
            return None, None
        # Fast sniff: extension or first line with tabs
        is_tsvish = tsv_path.suffix.lower() in {".tsv", ".txt"} or "\t" in tsv_path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
        if not is_tsvish:
            return None, None

        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        if not {"contig_id", "was_circular"}.issubset(df.columns):
            return None, None

        rows = df.loc[df["contig_id"] == target_contig]
        if rows.empty:
            # No exact target row → treat as not applicable
            return None, None

        was_circ_raw = rows["was_circular"].iloc[0]
        was_circ = coerce_bool(was_circ_raw)
        if was_circ is None:
            return None, None

        stats = "circular genome" if was_circ else "scaffold"
        avg_coverage = None
        if "avg_coverage" in rows.columns:
            cov_raw = rows["avg_coverage"].iloc[0]
            if cov_raw is not None and not (isinstance(cov_raw, float) and np.isnan(cov_raw)):
                cov_str = str(cov_raw).strip()
                if cov_str and cov_str.upper() != "NA":
                    avg_coverage = cov_str

        return stats, avg_coverage
    except Exception as e:
        print(f"⚠️ Could not parse contig_stats-like file: {e}")
        return None, None

def parse_log_for_stats_and_cov(log_text: str):
    """Parse GetOrganelle-style log for stats, avg_coverage, avg_base_coverage."""
    match_stats = re.findall(r"Result status of animal_mt:\s*(.+)", log_text)
    stats = match_stats[-1].strip() if match_stats else None

    # GetOrganelle labels this line "Average animal_mt coverage = ..." for
    # scaffold results but "Average animal_mt kmer-coverage = ..." for circular
    # genomes. Accept both so circular assemblies (e.g. successful reseeds) still
    # populate avg_coverage instead of recording NULL.
    match_avg_coverage = re.findall(r"Average animal_mt (?:kmer-)?coverage =\s*([^\s]+)", log_text)
    avg_coverage = match_avg_coverage[-1].strip() if match_avg_coverage else None

    match_avg_base_coverage = re.findall(r"Average animal_mt base-coverage =\s*([^\s]+)", log_text)
    avg_base_coverage = match_avg_base_coverage[-1].strip() if match_avg_base_coverage else None

    return stats, avg_coverage, avg_base_coverage

# -------------------------------
# Main logic
# -------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Push mitogenome assembly stats to the OceanOmics DB."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite an existing mitogenome_data row even when it already "
            "holds a real result. Off by default: existing successful rows "
            "are preserved and only NULL / 'failed to assemble' rows get "
            "updated automatically."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("assembly_prefix")
    parser.add_argument(
        "input_path",
        help="GetOrganelle log OR mitohifi contigs_stats.tsv.",
    )
    parser.add_argument("fasta_path")
    args = parser.parse_args()

    config_file = args.config_file
    assembly_prefix = args.assembly_prefix
    input_path = Path(args.input_path)
    fasta_path = Path(args.fasta_path)
    force_overwrite = args.force

    # Defaults
    stats = None
    avg_coverage = None
    avg_base_coverage = None
    failed_to_assemble = False

    # Compute sequence length from FASTA up-front so we can detect the
    # "process finished but produced no contig" case and short-circuit.
    try:
        with open(fasta_path) as f:
            length = sum(len(line.strip()) for line in f if not line.startswith(">"))
    except Exception as e:
        print(f"❌ Failed to read FASTA file: {e}")
        sys.exit(1)

    if length == 0:
        # Empty FASTA = assembler exited cleanly without producing a contig.
        failed_to_assemble = True
        stats = "failed to assemble"
        avg_coverage = None
        avg_base_coverage = None
        print("ℹ️ Empty assembly FASTA detected — recording 'failed to assemble'.")
    else:
        # First, try interpreting the 3rd arg as contig_stats.tsv
        stats_from_tsv, avg_coverage_from_tsv = try_parse_contig_stats(input_path, target_contig="final_mitogenome")

        if stats_from_tsv is not None:
            # HiFi-style input: use circular/scaffold mapping
            stats = stats_from_tsv
            avg_coverage = avg_coverage_from_tsv
            avg_base_coverage = avg_coverage_from_tsv
            print(f"ℹ️ Detected contig_stats.tsv. Using stats = '{stats}' from final_mitogenome row.")
        else:
            # Fall back to GetOrganelle log parsing
            try:
                log_text = input_path.read_text()
            except Exception as e:
                print(f"❌ Failed to read input file: {e}")
                sys.exit(1)

            stats, avg_coverage, avg_base_coverage = parse_log_for_stats_and_cov(log_text)
            print("ℹ️ Detected log file. Parsed GetOrganelle-style fields.")

    print(f"Stats: {stats}")
    print(f"Length: {length}")
    print(f"Avg Coverage: {avg_coverage}")
    print(f"Avg Base Coverage: {avg_base_coverage}")

    # Parse assembly_prefix
    try:
        og_id, tech, seq_date, code = assembly_prefix.split(".")
    except ValueError:
        print(f"❌ Failed to split assembly_prefix: {assembly_prefix}")
        sys.exit(1)

    try:
        db_params = load_db_config(config_file)
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()

        field_names = ["stats", "length", "avg_coverage", "avg_base_coverage"]
        params = {
            "og_id": og_id,
            "tech": tech,
            "seq_date": seq_date,
            "code": code,
            "stats": stats,
            "length": int(length),
            "avg_coverage": float(avg_coverage) if avg_coverage is not None else None,
            "avg_base_coverage": float(avg_base_coverage) if avg_base_coverage is not None else None,
        }

        # Look up the current row (if any) and decide whether to write.
        # Default policy is insert-only: existing rows are preserved unless
        #   - the new run is a successful assembly and the prior row holds
        #     NULL or 'failed to assemble' (a strict upgrade — never
        #     destructive), or
        #   - --force was passed.
        cursor.execute(
            """
            SELECT stats, length, avg_coverage, avg_base_coverage
            FROM mitogenome_data
            WHERE og_id = %(og_id)s
              AND tech = %(tech)s
              AND seq_date = %(seq_date)s
              AND code = %(code)s
            """,
            {"og_id": og_id, "tech": tech, "seq_date": seq_date, "code": code},
        )
        existing = cursor.fetchone()

        should_write = True
        skip_reason = None
        if existing is not None and not force_overwrite:
            existing_stats = existing[0]
            existing_is_failed = (
                existing_stats is not None
                and str(existing_stats).strip().lower() == "failed to assemble"
            )
            if existing_stats is None:
                # Incomplete prior row; overwriting is strictly an improvement.
                pass
            elif existing_is_failed and not failed_to_assemble:
                # New run upgrades a prior failure to a success — always allow.
                pass
            elif existing_is_failed and failed_to_assemble:
                should_write = False
                skip_reason = "row already records 'failed to assemble' (no change needed)"
            else:
                should_write = False
                if failed_to_assemble:
                    skip_reason = (
                        "existing successful row found — refusing to overwrite with "
                        "'failed to assemble'"
                    )
                else:
                    skip_reason = "row already exists; pass --force to overwrite"

        if not should_write:
            existing_dict = dict(zip(field_names, existing))
            print(f"⚠️ Existing values preserved for {assembly_prefix}: {skip_reason}.")
            print(f"📌 Preserved stored values: {existing_dict}")
        else:
            upsert_query = """
            INSERT INTO mitogenome_data (
                og_id, tech, seq_date, code, stats, length, avg_coverage, avg_base_coverage
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(stats)s, %(length)s, %(avg_coverage)s, %(avg_base_coverage)s
            )
            ON CONFLICT (og_id, tech, seq_date, code)
            DO UPDATE SET
                stats = EXCLUDED.stats,
                length = EXCLUDED.length,
                avg_coverage = EXCLUDED.avg_coverage,
                avg_base_coverage = EXCLUDED.avg_base_coverage
            RETURNING stats, length, avg_coverage, avg_base_coverage
            """
            cursor.execute(upsert_query, params)
            returned = cursor.fetchone()
            conn.commit()

            final_dict = dict(zip(field_names, returned))
            print(f"📌 Final stored values: {final_dict}")
            if force_overwrite and existing is not None:
                print(
                    f"✅ Success: Overwrote mitogenome_data for {assembly_prefix} (--force)."
                )
            else:
                print(f"✅ Success: Inserted/Updated mitogenome_data for {assembly_prefix}")

    except Exception as e:
        if 'conn' in locals():
            conn.rollback()
        print(f"❌ Database error: {e}")

    finally:
        if 'cursor' in locals():
            cursor.close()
        if 'conn' in locals():
            conn.close()
