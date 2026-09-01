#!/usr/bin/env python3
"""
Compile the per-sample SQL-upload status files into a single report.

Inputs are the ``*.upload.txt`` files produced by:

  * ``PUSH_MTDNA_ASSM_RESULTS``        -> ``<og_id>.mtdna.upload.txt``
  * ``PUSH_MTDNA_ANNOTATION_RESULTS``  -> ``<og_id>.annotation.upload.txt``
  * ``SPECIES_VALIDATION``             -> ``<og_id>.species_validation.upload.txt``
  * ``PUSH_LCA_BLAST_RESULTS``         -> ``<og_id>.lca_blast.upload.txt``
  * ``PUSH_ENA_VALIDATION_RESULTS``    -> ``<assembly_prefix>[.<annotation>].ena_validation.upload.txt``
  * ``PUSH_QC_VALIDATOR``              -> ``<assembly_prefix>[.<annotation>].qc_validator.upload.txt``

The last two are named after ``meta.full_seqid``, so they carry the annotation
version that the other five steps omit; ``resolve_ena_key`` strips it back off
so every half lands on the same report row.

The script writes two output files:

  * ``upload_results_summary.tsv``  -- one row per sample, one column per
    upload step, with a normalised status string.
  * ``upload_results_appendix.txt`` -- the original detailed text from each
    upload file, grouped by sample, so issues can be diagnosed without
    chasing per-step files in the publish dir.
"""

import argparse
import os
import re
import sys
from pathlib import Path


# Map the trailing portion of the filename (between the og_id and the
# ``.upload.txt`` tail) to the column name that appears in the summary TSV.
STEP_FROM_SUFFIX = {
    "mtdna": "assembly",
    "annotation": "annotation",
    "species_validation": "species_validation",
    "lca_blast": "lca_blast",
    "lca_raw": "lca_raw",
    "ena_validation": "ena_validation",
    "qc_validator": "qc_validator",
}

# Column order in the summary TSV.
STEP_ORDER = [
    "assembly", "annotation", "species_validation", "lca_blast", "lca_raw",
    "ena_validation", "qc_validator",
]


def parse_filename(path):
    """
    Extract (assembly_prefix, og_id, step) from an upload status filename.

    Files are named ``<assembly_prefix>.<step>.upload.txt`` where the assembly
    prefix is the unique per-assembly key (e.g. ``OG868.hic.250624.getorg1770``)
    and ``<step>`` is the trailing token (``mtdna``, ``lca_raw`` …). The step is
    therefore the *last* dot-delimited token before ``.upload.txt`` and the OG id
    is the *first* token. Older single-token names (``OG868.mtdna.upload.txt``)
    still parse correctly, with the prefix equal to the OG id.

    Returns ``(None, None, None)`` if the file doesn't match the expected pattern.
    """
    name = path.name
    if not name.endswith(".upload.txt"):
        return None, None, None
    stem = name[: -len(".upload.txt")]
    parts = stem.split(".")
    if len(parts) < 2:
        return None, None, None
    suffix = parts[-1]
    step = STEP_FROM_SUFFIX.get(suffix)
    if step is None:
        return None, None, None
    assembly_prefix = ".".join(parts[:-1])
    og_id = parts[0]
    return assembly_prefix, og_id, step


def classify_status(step, text):
    """
    Return a normalised status string for the upload step based on the
    captured stdout text of its push script.
    """
    if text is None:
        return "missing"
    body = text.strip()
    if not body:
        return "missing"

    has_db_error = "❌ Database error" in body
    has_failure_marker = "❌" in body

    if step == "assembly":
        # Specific preserved-row reasons take precedence over the generic
        # "preserved" status so the report carries more useful information.
        if has_db_error:
            return "failed"
        if "refusing to overwrite with 'failed to assemble'" in body:
            return "failed_kept_prior_success"
        if "row already records 'failed to assemble'" in body:
            return "failed_to_assemble"
        # A preserved row that nonetheless gained the uniform depth measurement.
        # Must be tested BEFORE the plain "preserved" match below, which would
        # otherwise swallow it and hide the fact that the row was updated at all.
        # The assembly stats themselves are still preserved in this case.
        if (
            "row already exists; pass --force to overwrite" in body
            and "Depth metrics added" in body
        ):
            return "preserved_depth_updated"
        if "row already exists; pass --force to overwrite" in body:
            return "preserved"
        if "failed to assemble" in body or "Empty assembly FASTA detected" in body:
            return "failed_to_assemble"
        if "Success: Overwrote mitogenome_data" in body:
            return "success_forced"
        if re.search(r"✅ Success: Inserted/Updated", body):
            return "success"
        if has_failure_marker:
            return "failed"
        return "unknown"

    if step == "annotation":
        if has_db_error:
            return "failed"
        if "row already exists; pass --force to overwrite" in body:
            return "preserved"
        if "Success: Overwrote mitogenome_data" in body:
            return "success_forced"
        if "⚠️ Existing values preserved" in body:
            return "preserved"
        if re.search(r"✅ Success: Inserted/Updated", body):
            return "success"
        if has_failure_marker:
            return "failed"
        return "unknown"

    if step == "species_validation":
        if "❌ Database error" in body:
            return "failed"
        if "has no nominal_species_id" in body:
            return "no_nominal_species"
        if "✅ Success: lca_validation overwritten (--force)" in body:
            return "validated_forced"
        if "✅ Success: lca_validation upserted" in body:
            return "validated"
        if "Existing values preserved for lca_validation" in body:
            return "preserved"
        if "Sample" in body and "not validated" in body:
            return "not_validated"
        if has_failure_marker:
            return "failed"
        return "unknown"

    if step == "lca_raw":
        # The script prints one summary line at the very bottom:
        #   "✅ Success: lca_raw_results upload complete for <sample> (<mode>) — ..."
        # or "⚠️ ... finished with errors ...".
        if "⚠️" in body and "finished with errors" in body:
            return "partial"
        if "✅ Success: lca_raw_results upload complete" in body and "force-overwrite" in body:
            return "success_forced"
        if "✅ Success: lca_raw_results upload complete" in body and "preserved" in body:
            return "preserved"
        if "✅ Success: lca_raw_results upload complete" in body:
            return "success"
        if has_db_error:
            return "failed"
        if has_failure_marker:
            return "failed"
        return "unknown"

    if step == "lca_blast":
        # The lca_blast script emits one line per upload phase:
        #   "✅ BLAST upload complete: X rows succeeded, Y failed"
        #   "✅ LCA upload complete:   X rows succeeded, Y failed"
        # Treat any non-zero failure count as a partial upload.
        failed_counts = re.findall(
            r"upload complete:\s*\d+\s*rows succeeded,\s*(\d+)\s*failed", body
        )
        if failed_counts and any(int(n) > 0 for n in failed_counts):
            return "partial"
        if "🎉 All processing complete." in body:
            return "success"
        if has_db_error:
            return "failed"
        if has_failure_marker:
            return "failed"
        return "unknown"

    if step == "ena_validation":
        exit_match = re.search(r"^UPLOAD_EXIT(?:=|\t)(\d+)\s*$", body, re.MULTILINE)
        if has_db_error or (exit_match and int(exit_match.group(1)) != 0):
            return "failed"
        if "✅ Success: inserted ENA validation attempt" in body:
            return "success"
        # push_ena_validation_results.py takes this branch whenever the upsert
        # hits an existing row, which is the normal case on a re-run.
        if "🔁 Updated ENA validation attempt" in body:
            return "success_updated"
        if "⚠️ Exact ENA validation result already recorded" in body:
            return "preserved"
        if has_failure_marker:
            return "failed"
        return "unknown"

    if step == "qc_validator":
        exit_match = re.search(r"^UPLOAD_EXIT(?:=|\t)(\d+)\s*$", body, re.MULTILINE)
        if has_db_error or (exit_match and int(exit_match.group(1)) != 0):
            return "failed"
        if "✅ Success: lca_validation validator_2 set to" in body:
            return "validator_2_set"
        # A validator_2 that is already filled in is the expected outcome for
        # any re-run, and for a sample a human signed off on first. It is never
        # overwritten, so this is a normal status and not a failure.
        if "⚠️ Existing validator_2 preserved" in body:
            return "preserved"
        if "not submission_ready" in body:
            return "not_ready"
        if "No lca_validation row" in body:
            return "no_row"
        if has_failure_marker:
            return "failed"
        return "unknown"

    return "unknown"


# Steps whose upload file is named after meta.full_seqid rather than
# meta.mt_assembly_prefix, so the annotation token has to be reconciled away
# before the row is grouped with the other steps'.
SEQID_KEYED_STEPS = {"ena_validation", "qc_validator"}

# Used only when two annotation versions of one assembly collapse onto a single
# row: the most severe status wins so a failure is never hidden by a sibling
# success. Anything not listed (success*, preserved*) ranks lowest.
STATUS_SEVERITY = {
    "failed": 5,
    "failed_to_assemble": 5,
    "failed_kept_prior_success": 5,
    "partial": 4,
    "unknown": 3,
    "missing": 2,
}


def _severity(status):
    return STATUS_SEVERITY.get(status, 1)


def _read_text(path):
    try:
        return path.read_text(errors="replace")
    except Exception as e:
        print(f"[WARN] Could not read {path}: {e}", file=sys.stderr)
        return None


def resolve_ena_key(raw_prefix, anchor_keys):
    """
    Map an ``ena_validation`` file's raw key onto the assembly prefix used by
    every other upload step, returning ``(assembly_prefix, annotation)``.

    PUSH_ENA_VALIDATION_RESULTS names its file after ``meta.full_seqid``, which
    is ``<mt_assembly_prefix>.<annotation_version>`` (for example
    ``OG868.hic.250624.v3mitohifi.emma102``), while the five other steps use the
    bare ``mt_assembly_prefix``. Without this reconciliation each assembly ends
    up as two half-empty rows in the summary table.

    ``anchor_keys`` is the set of prefixes seen from the other steps, so the
    annotation token is stripped only when doing so actually lands on a real
    assembly. When no other step ran (the standalone ``ena.nf`` entry point)
    fall back to the positional split used by ``identity_from_seqid`` in
    ``collate_ena_validation.py``: four fields of identity, the rest annotation.
    """
    if raw_prefix in anchor_keys:
        return raw_prefix, ""
    head, _, tail = raw_prefix.rpartition(".")
    if head and head in anchor_keys:
        return head, tail
    parts = raw_prefix.split(".")
    if len(parts) > 4:
        return ".".join(parts[:4]), ".".join(parts[4:])
    return raw_prefix, ""


def collect_inputs(input_dir):
    """
    Walk the input directory and group upload files by assembly prefix (the
    unique per-assembly key) so each assembly gets one report row even when a
    single OG has several assemblies. Returns
    {assembly_prefix: {"og_id": og_id, "annotation": version, "steps": {step: (path, text)}}}.
    """
    parsed = []
    for path in sorted(Path(input_dir).rglob("*.upload.txt")):
        raw_prefix, og_id, step = parse_filename(path)
        if raw_prefix is None or step is None:
            print(
                f"[WARN] Skipping unrecognised upload file: {path.name}",
                file=sys.stderr,
            )
            continue
        parsed.append((path, raw_prefix, og_id, step))

    # Prefixes contributed by the steps that already agree on mt_assembly_prefix.
    anchor_keys = {prefix for _, prefix, _, step in parsed if step not in SEQID_KEYED_STEPS}

    grouped = {}
    for path, raw_prefix, og_id, step in parsed:
        assembly_prefix, annotation = raw_prefix, ""
        if step in SEQID_KEYED_STEPS:
            assembly_prefix, annotation = resolve_ena_key(raw_prefix, anchor_keys)
        text = _read_text(path)
        entry = grouped.setdefault(
            assembly_prefix, {"og_id": og_id, "annotation": "", "steps": {}}
        )
        existing = entry["steps"].get(step)
        if existing is not None:
            # Two annotation versions of the same assembly. Keep the worse of
            # the two statuses so a failure is still visible in the table.
            print(
                f"[WARN] {assembly_prefix}: multiple {step} files collapsed onto "
                f"one row ({existing[0].name}, {path.name})",
                file=sys.stderr,
            )
            if _severity(classify_status(step, text)) > _severity(
                classify_status(step, existing[1])
            ):
                entry["steps"][step] = (path, text)
        else:
            entry["steps"][step] = (path, text)
        if annotation:
            versions = {v for v in entry["annotation"].split(",") if v}
            versions.add(annotation)
            entry["annotation"] = ",".join(sorted(versions))
    return grouped


def write_summary_tsv(grouped, output_path):
    header = ["assembly_prefix", "og_id", "annotation_version"] + STEP_ORDER
    with open(output_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for assembly_prefix in sorted(grouped):
            steps = grouped[assembly_prefix]["steps"]
            row = [
                assembly_prefix,
                grouped[assembly_prefix]["og_id"],
                grouped[assembly_prefix].get("annotation", ""),
            ]
            for step in STEP_ORDER:
                entry = steps.get(step)
                text = entry[1] if entry else None
                row.append(classify_status(step, text))
            f.write("\t".join(row) + "\n")


def write_appendix(grouped, output_path):
    with open(output_path, "w") as f:
        f.write("Per-sample SQL-upload report appendix\n")
        f.write("=" * 60 + "\n\n")
        for assembly_prefix in sorted(grouped):
            steps = grouped[assembly_prefix]["steps"]
            annotation = grouped[assembly_prefix].get("annotation", "")
            heading = assembly_prefix + (f" (annotation {annotation})" if annotation else "")
            f.write(f"## Assembly: {heading}\n")
            for step in STEP_ORDER:
                entry = steps.get(step)
                f.write(f"\n---- {assembly_prefix} :: {step} ----\n")
                if not entry:
                    f.write("(no upload file produced for this step)\n")
                    continue
                path, text = entry
                f.write(f"# source: {path.name}\n")
                if text is None:
                    f.write("(file unreadable)\n")
                else:
                    f.write(text.rstrip() + "\n")
            f.write("\n" + ("=" * 60) + "\n\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the per-sample <og_id>.<step>.upload.txt files.",
    )
    parser.add_argument(
        "--summary",
        default="upload_results_summary_mqc.tsv",
        help="Output path for the TSV summary table.",
    )
    parser.add_argument(
        "--appendix",
        default="upload_results_appendix.txt",
        help="Output path for the concatenated detail appendix.",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"❌ --input-dir not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    grouped = collect_inputs(args.input_dir)
    write_summary_tsv(grouped, args.summary)
    write_appendix(grouped, args.appendix)

    n_samples = len(grouped)
    print(f"✅ Wrote summary for {n_samples} sample(s) to {args.summary}")
    print(f"📝 Detailed appendix at {args.appendix}")


if __name__ == "__main__":
    main()
