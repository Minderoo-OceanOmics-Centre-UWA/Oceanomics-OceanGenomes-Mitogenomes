#!/usr/bin/env python3
"""
Snapshot a completed mitogenome run, so a later run can be compared to it
mechanically rather than from memory.

Written for the NOVA_260909_LA remediation: the point of the rerun is to show
that 22 ATP9 and 3 mtMutS annotations survive, that 15 of 32 "broken" CDS were
QC false positives, that three complete fish stop failing on gene order, and that
47 assemblies held with no stated reason acquire one. None of that is checkable
unless the "before" is recorded in the same shape as the "after".

Reads only published run outputs. Writes three TSVs plus a summary:

    <out>/assemblies.tsv   one row per assembly: QC verdict, annotation, holds
    <out>/cds_integrity.tsv one row per broken CDS, classified under the sample's
                            own genetic code
    <out>/raw_features.tsv  every raw MITOS feature and whether it survived
    <out>/summary.txt       the counts worth quoting

Usage:
    mito_run_snapshot.py --run <run dir> --out <snapshot dir>
"""

import argparse
import collections
import csv
import json
import glob
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from annotation_stats import score_cds_dir
from mitos_to_emma import EXTRA_PCG_MAP, NON_GENE_MAP

# Raw MITOS names that are not the vertebrate set, tracked explicitly because
# these are the ones the conversion used to drop.
LINEAGE_FEATURES = ("atp9", "msh1", "mutS", "mtMutS", "OH", "OL", "gpI", "gpII", "lagli")

# What each raw MITOS name is CALLED once mitos_to_emma.py has run. Searching the
# published GFF for the raw name reports every renamed feature as lost: msh1 is
# published as MTMUTS and OH as REP_ORIGIN_H, so the three mtMutS and 61 OH
# annotations this script once reported dropped were in fact all present. Only
# names the rename preserves case-insensitively (atp9 -> ATP9, lagli -> LAGLI)
# ever looked right. Both tables are imported rather than copied so a new gene
# added to the converter cannot silently regress this again.
PUBLISHED_NAME = {raw: internal for raw, internal in EXTRA_PCG_MAP.items()}
PUBLISHED_NAME.update({raw: spec[0] for raw, spec in NON_GENE_MAP.items()})


def published_name(bare_name):
    """The name `bare_name` is published under, after the conversion rename."""
    return PUBLISHED_NAME.get(bare_name.lower(), bare_name)


def samplesheet_from_params(run):
    """The --input the run was actually launched with, from pipeline_info.

    A run launched with --input has no <outdir>/samplesheet/ at all, because the
    sheet lives with whichever run created it. Without this the genetic codes
    resolve to nothing, every CDS is skipped, and the snapshot reports zero broken
    CDS -- which reads as "clean" and is indistinguishable from a real pass. That
    is the wrong direction for a tool whose whole job is to be compared against.
    """
    for params in sorted(Path(run).glob("pipeline_info/params_*.json"), reverse=True):
        try:
            with open(params) as fh:
                value = json.load(fh).get("input")
        except (OSError, ValueError):
            continue
        if value and Path(value).is_file():
            return Path(value)
    return None


def read_taxonomy(run, samplesheet=None):
    """{og_id: {class, genetic_code, ...}}, from whichever source the run has.

    taxonomy_resolution.tsv first, since it records how each class was resolved;
    the samplesheet actually in use as the fallback, which is the only source a
    --input run has. Returns {} when neither is available, and the caller reports
    the CDS scan as not performed rather than as clean.
    """
    path = Path(run) / "samplesheet" / "taxonomy_resolution.tsv"
    if path.is_file():
        with open(path) as fh:
            return {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}

    sheet = Path(samplesheet) if samplesheet else samplesheet_from_params(run)
    if not sheet or not sheet.is_file():
        return {}
    with open(sheet) as fh:
        rows = {}
        for r in csv.DictReader(fh):
            key = r.get("og_id") or r.get("sample")
            if key:
                rows[key] = r
        return rows


def read_fasta(path):
    seq, started = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def bare(name):
    return re.sub(r"[-_](?:\d+|[a-z])$", "", re.sub(r"\(.*?\)$", "", name.strip()))


def snapshot(run, out, samplesheet=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tax = read_taxonomy(run, samplesheet)

    # ---- QC verdicts -----------------------------------------------------
    assemblies = []
    for path in sorted(glob.glob(f"{run}/qc/*.qc_summary.tsv")):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                og = row.get("sample", "")
                t = tax.get(og, {})
                row["class"] = t.get("class", "")
                row["order"] = t.get("order", "")
                row["genetic_code"] = t.get("genetic_code", "")
                row["nominal_species_id"] = t.get("nominal_species_id", "")
                row["taxonomy_source"] = t.get("source", "")
                assemblies.append(row)

    # ---- CDS integrity ---------------------------------------------------
    # Scored by the pipeline's own scorer, not a local copy of the loop. The copy
    # had no terminator-convention detection, so every EMMA (fish) CDS was flagged
    # no_stop: 136 false findings across 9 samples, 13 of 13 PCGs each, because
    # EMMA publishes its per-gene CDS with the stop codon removed.
    broken = []
    for cds_dir in sorted(glob.glob(f"{run}/mitogenomes/*/*/annotation/cds")):
        parts = Path(cds_dir).parts
        og, assembly = parts[-4], parts[-3]
        code = (tax.get(og, {}) or {}).get("genetic_code") or ""
        if not code:
            continue
        code = int(code)
        lineage = (tax.get(og, {}) or {}).get("class", "")
        # Pass the GFF for the same reason the pipeline does: a stop-excluded annotator's
        # per-gene files cannot be terminator-checked, and the GFF span can. Without it this
        # snapshot would disagree with the gate it is meant to be compared against.
        gffs = sorted(glob.glob(f"{run}/mitogenomes/{og}/{assembly}/annotation/*.gff"))
        # The run-level convention note is deliberately discarded: this TSV is one row per
        # CDS, and an annotator-wide advisory is not a finding about any gene.
        per_gene, _stop_excluded, _note = score_cds_dir(cds_dir, code, lineage=lineage,
                                                       gff_path=gffs[0] if gffs else None)
        if not per_gene:
            continue
        for gene, gene_issues, gene_advisories in per_gene:
            # score_cds_dir prefixes each finding with the gene; this TSV has a
            # gene column already, so strip it back off.
            findings = [f.split(":", 1)[1] if f.startswith(f"{gene}:") else f
                        for f in gene_issues + gene_advisories]
            if findings:
                broken.append({
                    "og_id": og, "assembly": assembly, "gene": gene,
                    "genetic_code": code, "class": lineage,
                    "issues": ";".join(findings),
                })

    # ---- raw MITOS features and their fate -------------------------------
    features = []
    for bed in sorted(glob.glob(f"{run}/mitogenomes/*/*/annotation/mitos_raw/result.bed")):
        parts = Path(bed).parts
        og, assembly = parts[-5], parts[-4]
        gff = sorted(glob.glob(f"{run}/mitogenomes/{og}/{assembly}/annotation/*.gff"))
        published = open(gff[0]).read() if gff else ""
        for line in open(bed):
            cols = line.strip().split("\t")
            if len(cols) < 6:
                cols = line.strip().split()
            if len(cols) < 6:
                continue
            raw = cols[3]
            b = bare(raw)
            pub = published_name(b)
            survived = bool(re.search(rf"Name=(MT-)?{re.escape(pub)}[;\s]", published,
                                      re.I))
            features.append({
                "og_id": og, "assembly": assembly, "raw_name": raw,
                "bare_name": b, "published_name": pub,
                "is_lineage_feature": "yes" if b in LINEAGE_FEATURES else "no",
                "survived_to_published_gff": "yes" if survived else "no",
            })

    def write(name, rows, fields):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    if assemblies:
        write("assemblies.tsv", assemblies, list(assemblies[0].keys()))
    write("cds_integrity.tsv", broken,
          ["og_id", "assembly", "gene", "genetic_code", "class", "issues"])
    write("raw_features.tsv", features,
          ["og_id", "assembly", "raw_name", "bare_name", "published_name",
           "is_lineage_feature", "survived_to_published_gff"])

    # ---- summary ---------------------------------------------------------
    lost = collections.Counter()
    kept = collections.Counter()
    for f in features:
        if f["is_lineage_feature"] == "yes":
            (kept if f["survived_to_published_gff"] == "yes" else lost)[f["bare_name"]] += 1

    broken_samples = {(b["og_id"], b["assembly"]) for b in broken
                      if "alternative_start" not in b["issues"] or len(b["issues"].split(";")) > 1}
    held = [a for a in assemblies if a.get("proceed_qc") == "false"]
    held_clean = [a for a in held if a.get("annotation_passed") == "yes"]

    lines = [
        f"run: {run}",
        f"assemblies: {len(assemblies)}",
        f"proceed_qc=true: {sum(1 for a in assemblies if a.get('proceed_qc') == 'true')}",
        f"held: {len(held)} (of which annotation_passed=yes: {len(held_clean)})",
        f"assemblies with >=1 flagged CDS: {len(broken_samples)}",
        "  NOTE: CDS are classified with the CURRENT orf_utils tables, not the",
        "  tables in force when the run executed. So this count is 'still broken",
        "  under the corrected rules', which is the set the fixer has to address.",
        "  It is deliberately NOT the count the old code produced: that number was",
        "  inflated by the code-9 start table and the missing arthropod CO1",
        "  exception, and reproducing it would just re-record two known bugs.",
        f"raw MITOS features: {len(features)}",
        "",
        "lineage-specific features LOST in conversion:",
    ]
    lines += [f"  {k}: {v}" for k, v in sorted(lost.items())] or ["  (none)"]
    lines += ["", "lineage-specific features KEPT:"]
    lines += [f"  {k}: {v}" for k, v in sorted(kept.items())] or ["  (none)"]
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[mito_run_snapshot] wrote 3 TSVs + summary.txt to {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="pipeline --outdir of a completed run")
    ap.add_argument("--out", required=True, help="snapshot directory to write")
    ap.add_argument("--samplesheet", default=None,
                    help="samplesheet to take class/genetic_code from. Only needed "
                         "when the run has neither samplesheet/taxonomy_resolution.tsv "
                         "nor a readable --input in pipeline_info/params_*.json.")
    args = ap.parse_args()
    snapshot(args.run, args.out, args.samplesheet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
