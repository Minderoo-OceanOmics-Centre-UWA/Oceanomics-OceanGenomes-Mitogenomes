#!/usr/bin/env python3
"""Flag when the mitogenome reference chosen for a sample is not relevant to it.

The reference is resolved by MITOHIFI_FINDMITOREFERENCE from the sample's species
*label* (reference_species_id / nominal_species_id). A wrong or coarse label
therefore yields a wrong-family reference, which silently degrades both
GetOrganelle seeding and the anthozoan annotation fix (e.g. a Merulinid sample
labelled "Montipora grisea" gets an Acroporidae reference its 16S/nad5 cannot be
transferred from). MITOS/coral_fix then just emit a partial annotation with no
obvious cause.

This check is label-free and taxonomy-DB-free: it BLASTs the chosen reference
genome against the assembled mitogenome and measures how much of the assembly the
reference covers and at what identity. A relevant reference (same genus/family)
covers most of the molecule at high identity; a wrong-family reference aligns only
in short, low-identity patches. Observed separation on real coral data:

    same genus   ~100% covered, ~99.6% id
    same family   ~92% covered, ~96.9% id
    wrong family   ~13% covered, ~81.5% id   <- flagged

Emits one line:  PASS|MISMATCH|UNKNOWN \t <details>  and always exits 0, so a
missing/odd reference can never break the run -- it just records a review flag.

Usage:
    reference_relevance_check.py --assembly asm.fa --reference-gb ref.gb \
        --out PREFIX.reference_relevance.txt [--min-cov 0.60] [--min-pid 88.0]
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO


def ref_to_fasta(gb_path, fa_path):
    """Write the reference GenBank sequence(s) to FASTA; return (organism, family)."""
    recs = list(SeqIO.parse(str(gb_path), "genbank"))
    if not recs:
        raise ValueError("no records in reference GenBank")
    SeqIO.write(recs, str(fa_path), "fasta")
    rec = recs[0]
    organism = rec.annotations.get("organism", rec.id)
    lineage = rec.annotations.get("taxonomy", []) or []
    family = next((t for t in lineage if t.endswith("idae")), "")
    return organism, family


def assembly_length(fa_path):
    return sum(len(r.seq) for r in SeqIO.parse(str(fa_path), "fasta"))


def union_coverage(intervals):
    """Total length covered by a set of (lo, hi) inclusive intervals, no double count."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    covered = 0
    clo, chi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo > chi + 1:
            covered += chi - clo + 1
            clo, chi = lo, hi
        else:
            chi = max(chi, hi)
    covered += chi - clo + 1
    return covered


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True, type=Path)
    ap.add_argument("--reference-gb", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-cov", type=float, default=0.60,
                    help="Flag MISMATCH when the reference covers less than this "
                         "fraction of the assembly (default 0.60).")
    ap.add_argument("--min-pid", type=float, default=88.0,
                    help="Flag MISMATCH when the coverage-weighted mean identity is "
                         "below this (default 88.0).")
    args = ap.parse_args()

    def finish(state, msg):
        args.out.write_text(f"{state}\t{msg}\n")
        print(f"[reference_relevance] {state}: {msg}", file=sys.stderr)
        sys.exit(0)

    if not args.assembly.exists() or args.assembly.stat().st_size == 0:
        finish("UNKNOWN", f"assembly missing/empty: {args.assembly.name}")
    if not args.reference_gb.exists() or args.reference_gb.stat().st_size == 0:
        finish("UNKNOWN", f"reference missing/empty: {args.reference_gb.name}")

    try:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            ref_fa = fh.name
        organism, family = ref_to_fasta(args.reference_gb, ref_fa)
    except Exception as exc:
        finish("UNKNOWN", f"could not parse reference {args.reference_gb.name}: {exc}")

    asm_len = assembly_length(args.assembly)
    if asm_len == 0:
        finish("UNKNOWN", "assembly length 0")

    out = subprocess.run(
        ["blastn", "-query", ref_fa, "-subject", str(args.assembly),
         "-evalue", "1e-5", "-outfmt", "6 sstart send length pident"],
        capture_output=True, text=True).stdout
    Path(ref_fa).unlink(missing_ok=True)

    intervals, aln_bp, wpid = [], 0, 0.0
    for line in out.strip().splitlines():
        ss, se, ln, pid = line.split("\t")
        ss, se, ln, pid = int(ss), int(se), int(ln), float(pid)
        intervals.append((min(ss, se), max(ss, se)))
        aln_bp += ln
        wpid += ln * pid
    cov = union_coverage(intervals) / asm_len
    mean_pid = (wpid / aln_bp) if aln_bp else 0.0

    ref_desc = organism + (f" [{family}]" if family else "")
    detail = (f"cov={cov:.2f} pid={mean_pid:.1f} ref={ref_desc} "
              f"thresholds=cov>={args.min_cov:.2f},pid>={args.min_pid:.1f}")

    if aln_bp == 0:
        finish("MISMATCH", f"reference does not align to assembly; {detail}")
    if cov < args.min_cov or mean_pid < args.min_pid:
        finish("MISMATCH", detail)
    finish("PASS", detail)


if __name__ == "__main__":
    main()
