#!/usr/bin/env python3
"""Check that every samplesheet sample ended the run either held or submission-ready.

The hold accounting has three fragment sources, and each closes one known way a
sample can leave the pipeline. None of them closes the CLASS of defect, which is
that a sample can leave at a point nobody thought to instrument -- an ignored task
failure inside the assembly subworkflow drops a sample with no assembly in any
filtered channel, no QC row and no database row, and the run's exit status reports
it only as a two-line ignored-task warning naming the process, not the consequence.

This check needs to know nothing about why a sample is missing. It compares the
samplesheet against held + submission-ready and names whatever is in neither.

The expected set comes from the SAMPLESHEET, deliberately not from the database:
a sample lost before its mitogenome_data row was written is absent from the DB too,
so a DB-derived expected set would exclude it and this check would pass while the
sample was still missing.

Warn, do not abort. A reporting gap found at the end of a multi-day run should not
throw the run away; it should be impossible to miss in the log. Exit status is
always 0 for that reason.

Pure stdlib so it runs in the plain python container the module already uses.
"""

import argparse
import csv
import sys
from pathlib import Path


def read_samplesheet_ogs(path):
    """One OG id per line, as written by collectFile in the main workflow."""
    if not path or not Path(path).is_file():
        return []
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def read_held_ogs(path):
    """Sample ids from column 1 of held_samples.tsv, skipping its header.

    The file's grain is one row per assembly per stage, so a sample can appear more
    than once; a set is what this check wants.
    """
    ogs = set()
    if not path or not Path(path).is_file():
        return ogs
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return ogs
        for row in reader:
            if row and row[0].strip():
                ogs.add(row[0].strip())
    return ogs


def read_submission_ready_ogs(path):
    """og_id values from ena_run_summary.tsv whose submission_ready is true.

    Columns are resolved by NAME from the header rather than by position, so a
    column added upstream cannot silently shift what this reads. A file with no
    header, or without those columns, yields an empty set rather than raising --
    the check is a diagnostic and must not be the thing that fails a run.
    """
    ogs = set()
    if not path or not Path(path).is_file():
        return ogs
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "og_id" not in reader.fieldnames:
            return ogs
        has_flag = "submission_ready" in reader.fieldnames
        for row in reader:
            og = (row.get("og_id") or "").strip()
            if not og:
                continue
            if has_flag:
                if str(row.get("submission_ready") or "").strip().lower() == "true":
                    ogs.add(og)
            else:
                ogs.add(og)
    return ogs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet-ogs", required=True)
    parser.add_argument("--held", required=True)
    parser.add_argument("--ena-run-summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    expected = read_samplesheet_ogs(args.samplesheet_ogs)
    held = read_held_ogs(args.held)
    ready = read_submission_ready_ogs(args.ena_run_summary)

    # Order follows the samplesheet so the report is stable and reviewable.
    unaccounted = [og for og in expected if og not in held and og not in ready]

    lines = [
        f"samplesheet samples:     {len(set(expected))}",
        f"held (any stage):        {len(held)}",
        f"submission-ready:        {len(ready)}",
        f"unaccounted:             {len(unaccounted)}",
    ]

    if unaccounted:
        for og in unaccounted:
            warning = (
                f"WARNING: {og} is in the samplesheet but appears in neither "
                "held_samples.tsv nor the submission-ready set. A sample left the "
                "pipeline at a stage that emits no held fragment. This is a reporting "
                "gap: find where the sample left, either a filter with no fragment "
                "source or an ignored task failure upstream of one, and report it there."
            )
            print(warning, file=sys.stderr)
            lines.append(warning)
    else:
        lines.append(
            f"all {len(set(expected))} samplesheet sample(s) accounted for"
        )

    Path(args.output).write_text("\n".join(lines) + "\n")

    # Always 0. See the module header: this warns, it does not abort.
    return 0


if __name__ == "__main__":
    sys.exit(main())
