#!/usr/bin/env python3

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
    return "circular genome"/"scaffold" based on target row.
    If file missing or not a contig_stats file, return None.
    """
    try:
        if not tsv_path.exists():
            return None
        # Fast sniff: extension or first line with tabs
        is_tsvish = tsv_path.suffix.lower() in {".tsv", ".txt"} or "\t" in tsv_path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
        if not is_tsvish:
            return None

        df = pd.read_csv(tsv_path, sep="\t")
        if not {"contig_id", "was_circular"}.issubset(df.columns):
            return None

        rows = df.loc[df["contig_id"] == target_contig]
        if rows.empty:
            # No exact target row → treat as not applicable
            return None

        was_circ_raw = rows["was_circular"].iloc[0]
        was_circ = coerce_bool(was_circ_raw)
        if was_circ is None:
            return None

        return "circular genome" if was_circ else "scaffold"
    except Exception as e:
        print(f"⚠️ Could not parse contig_stats-like file: {e}")
        return None

def parse_log_for_stats_and_cov(log_text: str):
    """Parse GetOrganelle-style log for stats, avg_coverage, avg_base_coverage."""
    match_stats = re.findall(r"Result status of animal_mt:\s*(.+)", log_text)
    stats = match_stats[-1].strip() if match_stats else None

    match_avg_coverage = re.findall(r"Average animal_mt coverage =\s*([^\s]+)", log_text)
    avg_coverage = match_avg_coverage[-1].strip() if match_avg_coverage else None

    match_avg_base_coverage = re.findall(r"Average animal_mt base-coverage =\s*([^\s]+)", log_text)
    avg_base_coverage = match_avg_base_coverage[-1].strip() if match_avg_base_coverage else None

    return stats, avg_coverage, avg_base_coverage

# -------------------------------
# Main logic
# -------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage:\n  push_mtdna_assm_results.py <config_file> <assembly_prefix> <log_or_contigstats_file> <fasta_file>")
        sys.exit(1)

    config_file = sys.argv[1]
    assembly_prefix = sys.argv[2]
    input_path = Path(sys.argv[3])  # can be a log OR a contig_stats.tsv
    fasta_path = Path(sys.argv[4])

    # Defaults
    stats = None
    avg_coverage = None
    avg_base_coverage = None

    # First, try interpreting the 3rd arg as contig_stats.tsv
    stats_from_tsv = try_parse_contig_stats(input_path, target_contig="final_mitogenome")

    if stats_from_tsv is not None:
        # HiFi-style input: use circular/scaffold mapping
        stats = stats_from_tsv
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

    # Compute sequence length from FASTA
    try:
        with open(fasta_path) as f:
            length = sum(len(line.strip()) for line in f if not line.startswith(">"))
    except Exception as e:
        print(f"❌ Failed to read FASTA file: {e}")
        sys.exit(1)

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

        upsert_query = """
        INSERT INTO mitogenome_data (
            og_id, tech, seq_date, code, stats, length, avg_coverage, avg_base_coverage
        )
        VALUES (
            %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(stats)s, %(length)s, %(avg_coverage)s, %(avg_base_coverage)s
        )
        ON CONFLICT (og_id, tech, seq_date, code) 
        DO UPDATE SET
            stats = CASE WHEN mitogenome_data.stats IS NULL THEN EXCLUDED.stats ELSE mitogenome_data.stats END,
            length = CASE WHEN mitogenome_data.length IS NULL THEN EXCLUDED.length ELSE mitogenome_data.length END,
            avg_coverage = CASE WHEN mitogenome_data.avg_coverage IS NULL THEN EXCLUDED.avg_coverage ELSE mitogenome_data.avg_coverage END,
            avg_base_coverage = CASE WHEN mitogenome_data.avg_base_coverage IS NULL THEN EXCLUDED.avg_base_coverage ELSE mitogenome_data.avg_base_coverage END
        RETURNING stats, length, avg_coverage, avg_base_coverage
        """

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

        cursor.execute(upsert_query, params)
        returned = cursor.fetchone()
        conn.commit()

        print("📌 Final stored values:")
        print(dict(zip(["stats", "length", "avg_coverage", "avg_base_coverage"], returned)))

        # Field-by-field comparison
        field_names = ["stats", "length", "avg_coverage", "avg_base_coverage"]
        preserved_fields = []

        for idx, field in enumerate(field_names):
            if params[field] is None:
                continue  # Skipped intentionally
            elif returned[idx] != params[field]:
                preserved_fields.append(field)

        if preserved_fields:
            print(f"⚠️ Existing values preserved for: {', '.join(preserved_fields)}")
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
