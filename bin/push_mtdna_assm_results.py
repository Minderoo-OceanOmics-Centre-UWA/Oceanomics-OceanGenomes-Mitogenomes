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

# Columns written from the uniform remap-based depth measurement
# (MITOGENOME_COVERAGE / bin/mito_depth.py).
#
# This tuple is the ONLY thing the depth-only upgrade path below is allowed to
# write. It deliberately contains no assembly column: stats, length, avg_coverage
# and avg_base_coverage never appear here, in any code path, so adding a depth to
# a preserved row can never rewrite the assembly result it belongs to.
#
# avg_coverage / avg_base_coverage keep their historical (assembler-specific,
# mutually incomparable) meaning and are left exactly as they were. mean_depth is
# the number to use for any cross-platform comparison; depth_method says which of
# the two a given row carries.
DEPTH_COLUMNS = (
    "mean_depth",
    "median_depth",
    "depth_sd",
    "depth_cv",
    "breadth_1x",
    "breadth_10x",
    "mito_mapped_reads",
    "total_reads",
    "mito_read_fraction",
    "depth_target_length_bp",
    "depth_target_fasta",
    "depth_method",
)

# Written when this run produced no depth at all: the assembly never reached
# annotation (failed, under-length, or a discarded GetOrganelle variant), or the
# depth step was skipped. Distinct from the legacy_* labels the SQL migration
# stamps on rows that predate the uniform measurement entirely.
DEPTH_METHOD_NOT_MEASURED = "not_measured"


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

def try_parse_getorg_check(tsv_path: Path):
    """
    If tsv_path looks like a getorg_check.tsv circularity-evidence table -- the
    schema emitted by both OATK_CHECK (modules/local/oatk/check_circularity) and
    the GetOrganelle check, carrying a `final_verdict_circular` column -- map that
    corrected verdict to 'circular genome' / 'scaffold'.

    Returns None when the file is missing or is not a getorg_check-shaped TSV, so
    the caller can fall through to the free-text GetOrganelle log parser. Carries
    no coverage: this schema records topology only, so avg_coverage is left to the
    other parsers / the depth measurement.
    """
    try:
        if not tsv_path.exists():
            return None
        lines = tsv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
        if not lines or "\t" not in lines[0]:
            return None

        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        if df.empty or "final_verdict_circular" not in df.columns:
            return None

        verdict = coerce_bool(df["final_verdict_circular"].iloc[0])
        if verdict is None:
            return None
        return "circular genome" if verdict else "scaffold"
    except Exception as e:
        print(f"⚠️ Could not parse getorg_check-like file: {e}")
        return None

def _num(value):
    """TSV cell -> float, or None for blank / NA / unparseable."""
    if value is None:
        return None
    token = str(value).strip()
    if not token or token.upper() == "NA":
        return None
    try:
        return float(token)
    except ValueError:
        return None


def parse_depth_tsv(path):
    """Read the single data row of a <prefix>.mito_depth.tsv.

    Returns a dict keyed by DEPTH_COLUMNS, or None when there is no usable
    measurement -- which covers the header-only assets/placeholders/empty_mito_depth.tsv
    placeholder, a missing file, and a fail-open run of mito_depth.py. Callers
    treat None as "this row was never measured", never as "depth is zero".
    """
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        print(f"⚠️ Could not parse depth TSV {path}: {e}")
        return None
    if df.empty or "mean_depth" not in df.columns:
        return None

    row = df.iloc[0]
    mean_depth = _num(row.get("mean_depth"))
    if mean_depth is None:
        # A row with no mean depth carries no measurement worth recording.
        return None

    def text(column):
        if column not in df.columns:
            return None
        value = row.get(column)
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return None
        token = str(value).strip()
        return token if token and token.upper() != "NA" else None

    def integer(column):
        value = _num(row.get(column) if column in df.columns else None)
        return int(value) if value is not None else None

    return {
        "mean_depth": mean_depth,
        "median_depth": _num(row.get("median_depth")),
        "depth_sd": _num(row.get("sd_depth")),
        "depth_cv": _num(row.get("depth_cv")),
        "breadth_1x": _num(row.get("breadth_1x")),
        "breadth_10x": _num(row.get("breadth_10x")),
        "mito_mapped_reads": integer("mito_mapped_reads"),
        "total_reads": integer("total_reads"),
        "mito_read_fraction": _num(row.get("mito_read_fraction")),
        "depth_target_length_bp": integer("target_length_bp"),
        "depth_target_fasta": text("target_fasta"),
        "depth_method": text("depth_method") or "remap_full_v1",
    }


# Skip reasons emitted when an existing mitogenome_data row is preserved. Their
# exact text is parsed by compile_upload_report.py, so do not reword them without
# updating classify_status() there.
SKIP_ALREADY_FAILED = "row already records 'failed to assemble' (no change needed)"
SKIP_WOULD_DOWNGRADE = (
    "existing successful row found — refusing to overwrite with "
    "'failed to assemble'"
)
SKIP_ROW_EXISTS = "row already exists; pass --force to overwrite"


def decide_write(existing_stats, failed_to_assemble, force_overwrite):
    """Decide whether this run may write over an existing mitogenome_data row.

    Returns ``(should_write, skip_reason)``, with ``skip_reason`` None when the
    write goes ahead. ``existing_stats`` is the stored `stats` value, or None
    both for "no row at all" (handled by the caller) and for an incomplete row.

    Two of the three refusals hold regardless of ``force_overwrite``: replacing a
    successful assembly with 'failed to assemble' is data loss rather than a
    policy choice, and re-writing an already-failed row changes nothing. Only the
    plain "a real result is already here" case is what ``--force`` unlocks.
    """
    existing_is_failed = (
        existing_stats is not None
        and str(existing_stats).strip().lower() == "failed to assemble"
    )
    if existing_stats is None:
        # Incomplete prior row; overwriting is strictly an improvement.
        return True, None
    if existing_is_failed and not failed_to_assemble:
        # New run upgrades a prior failure to a success — always allow.
        return True, None
    if existing_is_failed:
        return False, SKIP_ALREADY_FAILED
    if failed_to_assemble:
        return False, SKIP_WOULD_DOWNGRADE
    if not force_overwrite:
        return False, SKIP_ROW_EXISTS
    return True, None


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
            "holds a successful result. Without it, such a row is preserved and "
            "only NULL / 'failed to assemble' rows get updated automatically. "
            "Never overrides the refusal to replace a successful row with "
            "'failed to assemble'."
        ),
    )
    parser.add_argument(
        "--depth-tsv",
        default=None,
        help=(
            "Uniform remap-based depth from MITOGENOME_COVERAGE "
            "(<prefix>.mito_depth.tsv). Optional: a header-only placeholder or an "
            "absent file records the row as depth_method='not_measured' rather "
            "than failing."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("assembly_prefix")
    parser.add_argument(
        "input_path",
        help="GetOrganelle log OR mitohifi contigs_stats.tsv.",
    )
    parser.add_argument("fasta_path")
    parser.add_argument(
        "--circular",
        default=None,
        help=(
            "Corrected circular verdict (true/false/null) from meta.circular, i.e. "
            "the GetOrganelle reference check / MitoHiFi circularity check. When a "
            "real (true/false) value is given it overrides the topology parsed from "
            "the log/contig_stats -- the check is authoritative for reseed/_rgj "
            "molecules GetOrganelle linearised and could not self-close. 'null' or "
            "unparseable leaves the log/tsv verdict untouched."
        ),
    )
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
    # Set when the assembly produced a contig but none of the parsers could
    # resolve a topology from the input. Carried down to the write so the run
    # stops claiming success on a row it left NULL (see below).
    stats_unresolved = False

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

        # Structured TSV first, free text last. Second attempt is the
        # getorg_check.tsv circularity-evidence schema that OATK_CHECK and the
        # GetOrganelle check both emit, so an oatk contig fed its evidence TSV
        # resolves a real topology instead of silently landing on stats=None.
        stats_from_check = try_parse_getorg_check(input_path)

        if stats_from_tsv is not None:
            # HiFi-style input: use circular/scaffold mapping
            stats = stats_from_tsv
            avg_coverage = avg_coverage_from_tsv
            avg_base_coverage = avg_coverage_from_tsv
            print(f"ℹ️ Detected contig_stats.tsv. Using stats = '{stats}' from final_mitogenome row.")
        elif stats_from_check is not None:
            # getorg_check.tsv evidence (oatk / GetOrganelle circularity check).
            stats = stats_from_check
            print(f"ℹ️ Detected getorg_check.tsv. Using stats = '{stats}' from final_verdict_circular.")
        else:
            # Fall back to GetOrganelle log parsing
            try:
                log_text = input_path.read_text()
            except Exception as e:
                print(f"❌ Failed to read input file: {e}")
                sys.exit(1)

            stats, avg_coverage, avg_base_coverage = parse_log_for_stats_and_cov(log_text)
            print("ℹ️ Detected log file. Parsed GetOrganelle-style fields.")

        # Apply the corrected circular verdict from meta.circular (the GetOrganelle
        # reference check). GetOrganelle linearises a closed molecule at an arbitrary
        # point and reports "N scaffold(s)"; the check confirms full reference
        # coverage and corrects the topology. This is the only signal carrying the
        # rgj (reference-guided-join) verdict, so honour it over the stale log.
        circular_override = coerce_bool(args.circular)
        if circular_override is not None:
            # A genuinely multi-scaffold result (GetOrganelle "N scaffold(s)" with
            # N > 1) is not a single closed molecule, so keep its descriptive
            # original stat rather than flattening it to a circular/scaffold label.
            # Only single-scaffold (or unspecified) results take the corrected
            # verdict -- this is the reseed/_rgj case the check actually closes.
            scaffold_match = re.search(r"(\d+)\s*scaffold", str(stats or ""), re.IGNORECASE)
            n_scaffold = int(scaffold_match.group(1)) if scaffold_match else None
            if n_scaffold is not None and n_scaffold > 1:
                print(f"ℹ️ Keeping original multi-scaffold stat '{stats}' (>1 scaffold); not applying --circular={args.circular}.")
            else:
                corrected = "circular genome" if circular_override else "scaffold"
                if corrected != stats:
                    print(f"ℹ️ Overriding stats '{stats}' -> '{corrected}' from --circular={args.circular}.")
                stats = corrected

        # An assembled contig (length > 0) whose topology none of the parsers and
        # no --circular override could resolve. The row is still written so the run
        # does not abort, but flag it loudly here and downgrade the success message
        # at the write so the NULL is greppable instead of invisible.
        if stats is None:
            stats_unresolved = True
            print(f"⚠️ Unrecognised assembly input: {input_path}")
            print("   Not a contigs_stats.tsv, getorg_check.tsv, or GetOrganelle log.")
            print("   Writing row with stats=NULL — needs manual review.")

    # Uniform remap-based depth. None when this assembly never reached annotation
    # (failed / under-length / a discarded variant) or the depth step was skipped.
    depth = parse_depth_tsv(args.depth_tsv)

    print(f"Stats: {stats}")
    print(f"Length: {length}")
    print(f"Avg Coverage: {avg_coverage}  (legacy, assembler-specific)")
    print(f"Avg Base Coverage: {avg_base_coverage}  (legacy, assembler-specific)")
    if depth is not None:
        print(f"Mean Depth: {depth['mean_depth']}  (uniform, {depth['depth_method']})")
    else:
        print("Mean Depth: None (not measured)")

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
        # Every row gets a depth_method so a legacy number is never mistaken for a
        # measured one, even when this run had nothing to measure.
        depth_values = dict(depth) if depth is not None else {
            column: None for column in DEPTH_COLUMNS
        }
        if depth is None:
            depth_values["depth_method"] = DEPTH_METHOD_NOT_MEASURED
        params.update(depth_values)

        # Look up the current row (if any) and decide whether to write.
        # Default policy is overwrite: an existing row is replaced by this run's
        # result. Two refusals hold regardless of --force, because neither is a
        # policy choice, they are the cases where writing would destroy data:
        #   - the prior row is a successful assembly and this run failed to
        #     assemble (never downgrade a result to 'failed to assemble');
        #   - the prior row already records 'failed to assemble' and so does this
        #     run (nothing to change).
        # Without --force one further case is preserved: a prior successful row
        # is kept when this run also succeeded, which is the insert-only
        # behaviour. A prior row with NULL stats is always overwritten, since
        # completing it is strictly an improvement.
        #
        # Phase R (docs/assembly_relations_plan.md §3): the assembly result lives in
        # mito_assembly and the depth verdict in mito_assembly_depth, each keyed on
        # assembly_id. Table names are unqualified on purpose: before the R.6 merge a test
        # session reaches the v3 tables with search_path = v3, public; after it they are in
        # public. The policy below is unchanged; only where each value lands has changed.
        cursor.execute(
            """
            SELECT stats, length, avg_coverage, avg_base_coverage, assembly_id
            FROM mito_assembly
            WHERE og_id = %(og_id)s
              AND tech = %(tech)s
              AND seq_date = %(seq_date)s
              AND code = %(code)s
            """,
            {"og_id": og_id, "tech": tech, "seq_date": seq_date, "code": code},
        )
        existing_row = cursor.fetchone()
        existing = existing_row[:4] if existing_row is not None else None

        should_write = True
        skip_reason = None
        if existing is not None:
            should_write, skip_reason = decide_write(
                existing[0], failed_to_assemble, force_overwrite
            )

        depth_columns_sql = ", ".join(DEPTH_COLUMNS)
        depth_values_sql = ", ".join("%({0})s".format(c) for c in DEPTH_COLUMNS)

        if not should_write:
            existing_dict = dict(zip(field_names, existing))
            print(f"⚠️ Existing values preserved for {assembly_prefix}: {skip_reason}.")
            print(f"📌 Preserved stored values: {existing_dict}")

            # Narrow, additive exception to the preserve-everything rule: an assembly
            # written before the uniform depth measurement existed has a perfectly good
            # assembly result but no mean_depth, and the guard above keys only on
            # `stats`, so it would otherwise never gain one without --force.
            #
            # This can only ever ADD a depth to an assembly that has none:
            #   1. it writes mito_assembly_depth only, so stats/length/avg_coverage/
            #      avg_base_coverage (in mito_assembly) cannot be touched from here;
            #   2. an existing depth row is replaced only WHERE its mean_depth IS NULL,
            #      so an assembly that already carries a measured depth is untouchable
            #      from here. Not ORed with the depth_method states: a row with a real
            #      mean_depth but a NULL / stale label would then have been eligible
            #      and its measurement overwritten. not_measured rows have mean_depth
            #      NULL, so this single predicate already covers them;
            #   3. it runs only when this run HAS a depth, so it can never null one out.
            if depth is not None:
                assignments = ", ".join(
                    "{0} = EXCLUDED.{0}".format(column) for column in DEPTH_COLUMNS
                )
                cursor.execute(
                    f"""
                    INSERT INTO mito_assembly_depth (assembly_id, {depth_columns_sql}, depth_measured_at)
                    VALUES (%(assembly_id)s, {depth_values_sql}, now())
                    ON CONFLICT (assembly_id) DO UPDATE
                    SET {assignments}, depth_measured_at = now()
                    WHERE mito_assembly_depth.mean_depth IS NULL
                    RETURNING mean_depth, depth_cv, depth_method
                    """,
                    {**params, "assembly_id": existing_row[4]},
                )
                depth_returned = cursor.fetchone()
                conn.commit()
                if depth_returned:
                    print(
                        f"📌 Depth metrics added for {assembly_prefix}: "
                        f"mean_depth={depth_returned[0]}, depth_cv={depth_returned[1]}, "
                        f"depth_method={depth_returned[2]}"
                    )
                    print(
                        f"✅ Depth metrics added for {assembly_prefix} "
                        "(assembly stats preserved)."
                    )
                else:
                    print(
                        f"ℹ️ Row for {assembly_prefix} already carries a measured depth; "
                        "left untouched."
                    )
        else:
            cursor.execute(
                """
                -- og_num is not stored on mito_assembly; the compatibility view derives it.
                INSERT INTO mito_assembly (
                    og_id, tech, seq_date, code, stats, length, avg_coverage, avg_base_coverage
                )
                VALUES (
                    %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(stats)s, %(length)s,
                    %(avg_coverage)s, %(avg_base_coverage)s
                )
                ON CONFLICT (og_id, tech, seq_date, code)
                DO UPDATE SET
                    stats = EXCLUDED.stats,
                    length = EXCLUDED.length,
                    avg_coverage = EXCLUDED.avg_coverage,
                    avg_base_coverage = EXCLUDED.avg_base_coverage,
                    updated_at = now()
                RETURNING assembly_id, stats, length, avg_coverage, avg_base_coverage
                """,
                params,
            )
            returned = cursor.fetchone()
            assembly_id = returned[0]

            # Replace the depth ONLY when this run actually measured one. Without this
            # guard a --force rerun that skipped the depth step (or whose assembly never
            # reached annotation) would null out a perfectly good existing measurement,
            # which is data loss rather than an overwrite. A new assembly with nothing
            # measured records not_measured, so the depth step's verdict is never
            # missing for an assembly this script wrote.
            if depth is not None:
                assignments = ", ".join(
                    "{0} = EXCLUDED.{0}".format(column) for column in DEPTH_COLUMNS
                )
                cursor.execute(
                    f"""
                    INSERT INTO mito_assembly_depth (assembly_id, {depth_columns_sql}, depth_measured_at)
                    VALUES (%(assembly_id)s, {depth_values_sql}, now())
                    ON CONFLICT (assembly_id) DO UPDATE
                    SET {assignments}, depth_measured_at = now()
                    """,
                    {**params, "assembly_id": assembly_id},
                )
            elif existing is None:
                cursor.execute(
                    """
                    INSERT INTO mito_assembly_depth (assembly_id, depth_method, depth_measured_at)
                    VALUES (%(assembly_id)s, %(depth_method)s, now())
                    ON CONFLICT (assembly_id) DO NOTHING
                    """,
                    {"assembly_id": assembly_id, "depth_method": DEPTH_METHOD_NOT_MEASURED},
                )
            cursor.execute(
                "SELECT mean_depth, depth_method FROM mito_assembly_depth WHERE assembly_id = %s",
                (assembly_id,),
            )
            depth_row = cursor.fetchone() or (None, None)
            conn.commit()

            final_dict = dict(zip(field_names + ["mean_depth", "depth_method"],
                                  tuple(returned[1:]) + tuple(depth_row)))
            print(f"📌 Final stored values: {final_dict}")
            # The status lines below are parsed by bin/compile_upload_report.py, which
            # is why they still say mitogenome_data. Change them together or not at all.
            if stats_unresolved:
                # A contig was assembled but its topology never resolved, so this
                # is not a success no matter how the write went. Say so, and make
                # it greppable across sql_uploaded_data/*.mtdna.upload.txt.
                print(f"⚠️ Wrote INCOMPLETE row for {assembly_prefix} (stats=NULL) — needs manual review.")
            elif force_overwrite and existing is not None:
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
