#!/usr/bin/env python3
"""Pick the per-sample coral annotation reference from the curated DB by SEQUENCE.

Replaces the label-based reference lookup (findMitoReference by species name) for
the coral/invert branch: BLAST the assembly against the bundled Anthozoa
mitogenome DB and emit the best-matching record's GenBank as the reference for
coral_fix_bed.py. This is label-free, so a wrong/coarse species label can no
longer hand a wrong-family reference to the annotation fixer.

Selection = the DB record with the highest total bitscore against the assembly
(closest + most-covering). The chosen record is reported with its assembly
coverage and coverage-weighted identity; below --min-cov/--min-pid it is still
emitted but marked LOW_CONFIDENCE so the downstream REFERENCE_RELEVANCE check /
QC flag can surface it (e.g. a lineage with no DB representative).

Always exits 0. Runs in the MITOS2 BioContainer (blastn + biopython).

Usage:
    select_coral_reference.py --assembly asm.fa --db-gb coral_mito_refdb.gb \
        --out-gb PREFIX.reference.gb --out-status PREFIX.reference_select.txt
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO


def union_coverage(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    covered, clo, chi = 0, *intervals[0]
    for lo, hi in intervals[1:]:
        if lo > chi + 1:
            covered += chi - clo + 1
            clo, chi = lo, hi
        else:
            chi = max(chi, hi)
    return covered + (chi - clo + 1)


def family_of(rec):
    return next((t for t in rec.annotations.get("taxonomy", []) or [] if t.endswith("idae")), "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True, type=Path)
    ap.add_argument("--db-gb", required=True, type=Path, help="coral_mito_refdb.gb (multi-record)")
    ap.add_argument("--out-gb", required=True, type=Path)
    ap.add_argument("--out-status", required=True, type=Path)
    ap.add_argument("--min-cov", type=float, default=0.60)
    ap.add_argument("--min-pid", type=float, default=88.0)
    args = ap.parse_args()

    def status(state, msg):
        args.out_status.write_text(f"{state}\t{msg}\n")
        print(f"[select_coral_reference] {state}: {msg}", file=sys.stderr)
        sys.exit(0)

    asm_recs = list(SeqIO.parse(str(args.assembly), "fasta"))
    if not asm_recs:
        status("NONE", f"empty assembly {args.assembly.name}")
    asm_len = sum(len(r.seq) for r in asm_recs)

    db = {r.id: r for r in SeqIO.parse(str(args.db_gb), "genbank")}
    if not db:
        status("NONE", f"empty DB {args.db_gb.name}")

    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        SeqIO.write(db.values(), fh.name, "fasta")
        db_fa = fh.name

    out = subprocess.run(
        ["blastn", "-query", str(args.assembly), "-subject", db_fa,
         "-evalue", "1e-10", "-max_hsps", "20",
         "-outfmt", "6 sseqid qstart qend pident length bitscore"],
        capture_output=True, text=True).stdout
    Path(db_fa).unlink(missing_ok=True)

    # Aggregate per DB subject: total bitscore (selection), assembly coverage +
    # weighted identity (the guard on the winner).
    agg = {}
    for line in out.strip().splitlines():
        sid, qs, qe, pid, ln, bits = line.split("\t")
        qs, qe, pid, ln, bits = int(qs), int(qe), float(pid), int(ln), float(bits)
        a = agg.setdefault(sid, {"bits": 0.0, "iv": [], "alnbp": 0, "wpid": 0.0})
        a["bits"] += bits
        a["iv"].append((min(qs, qe), max(qs, qe)))
        a["alnbp"] += ln
        a["wpid"] += ln * pid
    if not agg:
        status("NONE", "no DB record aligned to the assembly")

    best_id = max(agg, key=lambda k: agg[k]["bits"])
    a = agg[best_id]
    cov = union_coverage(a["iv"]) / asm_len
    pid = a["wpid"] / a["alnbp"] if a["alnbp"] else 0.0

    rec = db[best_id]
    # Re-id so coral_fix's status/provenance shows the accession, keep features.
    SeqIO.write([rec], str(args.out_gb), "genbank")

    organism = rec.annotations.get("organism", best_id)
    fam = family_of(rec)
    confident = cov >= args.min_cov and pid >= args.min_pid
    state = "SELECTED" if confident else "SELECTED_LOW_CONFIDENCE"
    status(state, f"{best_id} {organism} [{fam}] cov={cov:.2f} pid={pid:.1f}")


if __name__ == "__main__":
    main()
