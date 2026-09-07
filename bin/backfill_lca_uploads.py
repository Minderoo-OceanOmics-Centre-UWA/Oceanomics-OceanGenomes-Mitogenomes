#!/usr/bin/env python3
"""Re-upload LCA results that a finished run published but never got into the DB.

Everything the upload needs survives in the run's published output, so a lost
upload does not need the pipeline re-run: this rebuilds the two combined inputs
that SPECIES_VALIDATION would have produced and hands them to the same
``push_lca_blast_results.py`` / ``push_lca_raw_results.py`` the pipeline uses.

Inputs come from ``mitogenomes/<OG>/<assembly>/lca/``, not ``species_validation/``,
because an archived run may have had its ``work``, ``qc`` and
``species_validation`` directories pruned while the per-region LCA outputs
remain. The per-region ``blast.*.filtered.tsv`` carry the same 24 columns, in the
same order, as the ``blast_combined.*.tsv`` the pipeline builds from them.

Which push runs for an assembly is decided by audit_lca_db_coverage.py: only the
tables that are actually empty are pushed, so re-running this is safe once the
gaps are filled -- it then finds nothing to do.

``--force`` is deliberately not passed through. The rows being restored are
missing, so there is no superseded history to prune, and skipping the prune
removes any chance of deleting a good older result.

Usage:
    backfill_lca_uploads.py --log-dir <dir> <config.cfg> <run_dir> [<run_dir> ...]
    backfill_lca_uploads.py --dry-run <config.cfg> /scratch/.../batch-*
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from audit_lca_db_coverage import find_gaps, load_db_config, load_db_counts

# Imported rather than reimplemented: the combined files have to be byte-for-byte
# what SPECIES_VALIDATION would have written, or the push scripts parse them
# differently.
from species_validation import concatenate_files, concatenate_lca_files

BIN_DIR = Path(__file__).resolve().parent


def build_inputs(assembly_dir, workdir):
    """Rebuild lca_combined / blast_combined for one assembly. Returns their paths."""
    lca_dir = assembly_dir / "lca"
    prefix = assembly_dir.name
    lca_combined = workdir / f"lca_combined.{prefix}.tsv"
    blast_combined = workdir / f"blast_combined.{prefix}.tsv"
    concatenate_lca_files(sorted(lca_dir.glob("lca.*.tsv")), lca_combined)
    concatenate_files(sorted(lca_dir.glob("blast.*.filtered.tsv")), blast_combined)
    return lca_combined, blast_combined


def run_push(command, log_path):
    """Run one push script, tee-ing its output to a log the way the pipeline does."""
    print("   $ " + " ".join(str(c) for c in command))
    result = subprocess.run(
        [str(c) for c in command], capture_output=True, text=True
    )
    output = result.stdout + result.stderr
    if log_path:
        with open(log_path, "a") as handle:
            handle.write("$ " + " ".join(str(c) for c in command) + "\n")
            handle.write(output)
    for line in output.splitlines():
        print(f"     {line}")
    return result.returncode


def backfill_one(gap, config_file, log_dir, dry_run):
    assembly_dir = Path(gap["assembly_dir"])
    og_id = gap["og_id"]
    missing = gap["missing"]
    log_path = Path(log_dir) / f"{assembly_dir.name}.backfill.txt" if log_dir else None

    print(f"-- {gap['run']}/{assembly_dir.name}: missing {', '.join(sorted(missing))}")
    if dry_run:
        return 0

    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        # push_lca_blast_results.py uploads blast_filtered_lca and lca together
        # and is the only way to write either, so it runs if either is missing.
        if "blast_filtered_lca" in missing or "lca" in missing:
            lca_combined, blast_combined = build_inputs(assembly_dir, workdir)
            failures += run_push(
                [
                    BIN_DIR / "push_lca_blast_results.py",
                    config_file,
                    og_id,
                    lca_combined,
                    blast_combined,
                ],
                log_path,
            )
        if "lca_raw_results" in missing:
            failures += run_push(
                [BIN_DIR / "push_lca_raw_results.py", config_file, og_id]
                + sorted((assembly_dir / "lca").glob("lca_raw.*.tsv")),
                log_path,
            )
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file")
    parser.add_argument("run_dirs", nargs="+", help="Pipeline outdirs to backfill.")
    parser.add_argument("--log-dir", help="Write one <assembly>.backfill.txt per assembly here.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List what would be pushed without touching the database.",
    )
    args = parser.parse_args()

    if args.log_dir:
        Path(args.log_dir).mkdir(parents=True, exist_ok=True)

    gaps = find_gaps(args.run_dirs, load_db_counts(load_db_config(args.config_file)))
    if not gaps:
        print("Nothing to backfill: every published LCA result is in the database.")
        return 0

    print(f"{len(gaps)} assemblies to backfill\n")
    failed = [
        gap
        for gap in gaps
        if backfill_one(gap, args.config_file, args.log_dir, args.dry_run)
    ]

    print(f"\nbackfilled={len(gaps) - len(failed)} failed={len(failed)}")
    for gap in failed:
        print(f"   FAILED {gap['run']}/{gap['assembly']}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
