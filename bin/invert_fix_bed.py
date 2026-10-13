#!/usr/bin/env python3
"""
Reference-free ORF boundary repair for any invertebrate MITOS annotation.

The coral fixer (bin/coral_fix_bed.py) does two separable things, and only one of
them is about corals: it transfers rrnL/nad5/cox1/mtMutS from an Anthozoa GenBank
reference, and it snaps every PCG's boundaries to the nearest clean ORF using the
assembly alone. The second needs no reference, no taxonomy and no coral, but it
was reachable only down the coral branch.

In NOVA_260909_LA that cost 17 assemblies. Every one of them had a genuinely
broken CDS; not one was a coral, so not one could reach the repair. 13 were
Porifera, and five Demospongiae shared the identical defect (ND4L opening on TTT
and running to a sense codon) -- one systematic MITOS boundary miscall rather
than five independent bad assemblies.

This script wraps bin/orf_snap.py and nothing else. It never consults a
reference, so it cannot import a gene the assembly does not carry: the worst it
can do is decline to move a boundary. Where it declines, it says why, and the
downstream QC step holds the sample rather than publishing a broken CDS.

Usage:
    invert_fix_bed.py --bed result.bed --genome rotated.fa --code 5 \
        --out-bed result.fixed.bed --status invert_fix_status.tsv \
        [--lineage Malacostraca]
"""

import argparse
import sys
from pathlib import Path

from orf_snap import bare_name, parse_bed, snap_pcg_rows


def read_genome(path):
    """First FASTA record's sequence, uppercased. Stdlib only: this runs in the
    MITOS container alongside coral_fix_bed.py."""
    seq = []
    started = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def write_bed(path, rows):
    with open(path, "w") as out:
        for r in rows:
            out.write("\t".join(str(c) for c in r) + "\n")


def write_status(path, status, actions, fixed, unrepaired):
    """One status line plus a per-gene TSV of what moved and what did not.

    Mirrors the coral fixer's status artifact so the two are read the same way.
    Every rejected move is recorded with its reason: a boundary this script
    declined to touch is a fact the QC step needs, not an absence.
    """
    with open(path, "w") as out:
        out.write(f"{status}\t{'; '.join(actions) if actions else 'no changes'}\n")
        out.write("#gene\toutcome\tdetail\n")
        for gene, detail in fixed:
            out.write(f"{gene}\tfixed\t{detail}\n")
        for gene, detail in unrepaired:
            out.write(f"{gene}\tunrepaired\t{detail}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, type=Path, help="MITOS2 result.bed")
    ap.add_argument("--genome", required=True, type=Path, help="Rotated MITOS input FASTA")
    ap.add_argument("--out-bed", required=True, type=Path)
    ap.add_argument("--status", required=True, type=Path)
    ap.add_argument("--code", type=int, required=True,
                    help="Mito genetic code from meta.genetic_code. No default: a "
                         "wrong table silently redefines every start and stop.")
    ap.add_argument("--lineage", default=None,
                    help="Taxonomic class or phylum. Only widens the accepted "
                         "initiation codons (orf_utils.LINEAGE_ALT_STARTS), so "
                         "omitting it can make this script attempt a repair it "
                         "need not, never the reverse.")
    ap.add_argument("--window", type=int, default=36,
                    help="Max nt a boundary may move in either direction.")
    ap.add_argument("--len-tol", type=float, default=0.15,
                    help="Max fractional length change a snapped CDS may take.")
    ap.add_argument("--linear", action="store_true",
                    help="Molecule is not circular: never snap a boundary across its ends.")
    args = ap.parse_args()

    if not args.bed.exists():
        sys.exit(f"BED not found: {args.bed}")
    if not args.genome.exists():
        sys.exit(f"Genome not found: {args.genome}")

    rows = parse_bed(args.bed)
    gseq = read_genome(args.genome)
    if not rows or not gseq:
        write_bed(args.out_bed, rows)
        write_status(args.status, "SKIPPED",
                     ["empty BED" if not rows else "empty genome"], [], [])
        return 0

    fixes, actions, unrepaired = snap_pcg_rows(
        rows, gseq, args.code, window=args.window, len_tol=args.len_tol,
        lineage=args.lineage, circular=not args.linear)

    fixed_detail = []
    for i, (new_s0, new_e1) in fixes.items():
        gene = bare_name(rows[i][3])
        fixed_detail.append((gene, f"{rows[i][1]}-{rows[i][2]} -> {new_s0}-{new_e1}"))
        rows[i][1], rows[i][2] = str(new_s0), str(new_e1)

    unrepaired_detail = []
    for gene in unrepaired:
        why = next((a for a in actions if a.startswith(f"PCG-skip {gene} ")), "")
        unrepaired_detail.append((gene, why.split("(", 1)[-1].rstrip(")") if why else "unrepaired"))

    # PARTIAL, not FIXED, whenever anything was left broken: a sample with one
    # repaired gene and one still-broken gene has not been fixed, and calling it
    # FIXED would let it through the gate on the strength of the half that worked.
    if unrepaired:
        status = "PARTIAL" if fixes else "FAILED"
    else:
        status = "FIXED" if fixes else "SKIPPED"

    write_bed(args.out_bed, rows)
    write_status(args.status, status, actions, fixed_detail, unrepaired_detail)
    print(f"[invert_fix_bed] {status}: {len(fixes)} boundary fix(es), "
          f"{len(unrepaired)} left unrepaired")
    # Always exit 0: this is a repair attempt, not a gate. The QC step decides
    # whether what came out is publishable.
    return 0


if __name__ == "__main__":
    sys.exit(main())
