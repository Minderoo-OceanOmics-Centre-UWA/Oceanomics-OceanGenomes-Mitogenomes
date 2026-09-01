#!/usr/bin/env python3
"""Re-origin a circular mitogenome assembly to the start of cox1 (intermediate).

Note: cox1 is an *intermediate* origin, not the final one. This rotation exists
only to keep MITOS2 from splitting the coral nad5 group I intron across the
linearisation point (see below). After MITOS2 annotates, mitos_to_emma.py
re-origins the published genome again -- to tRNA-Met's 5' end -- so the deposited
mitogenome starts at trnM, the convention NCBI coral genomes follow. trnM is
located from MITOS's own annotation (no BLAST/tblastn: a ~70 bp tRNA is conserved
in structure, not enough in primary sequence, for a reliable nucleotide search).

Why: MITOS2 annotates the invertebrate (coral) path. Scleractinian/anthozoan
nad5 carries a large group I intron that, with most other genes nested inside
it, wraps GetOrganelle's arbitrary linearisation point. When that wrap crosses
position 1, MITOS splits nad5 across the origin and emits coordinates whose
linear span swallows unrelated genes (cox1, rrnL, ...), which then looks like a
"missing nad5" downstream and produces overlapping features that NCBI/table2asn
reject. Rotating the assembly so position 1 is the start of cox1 -- a conserved,
single-exon gene corals always carry -- puts every feature in a sane linear
order with nothing crossing the origin. This mirrors EMMA's `--rotate MT-TF` for
vertebrates (corals lack tRNA-Phe, so cox1 is the anchor instead).

cox1 is located with tblastn against a small panel of coral cox1 proteins. The
rotation only needs to land near cox1's 5' end, and a few bp of imprecision is
harmless.

cox1 is usually single-exon, but not always: some scleractinians carry a second
group I intron in cox1 (holding a LAGLIDADG homing endonuclease ORF), so cox1
comes back as two tblastn HSPs on the same diagonal, split by ~1 kb of intron.
Anchoring on the best-scoring HSP is wrong there -- the 3' exon usually scores
higher, and back-extrapolating 3*(qstart-1) from it lands the origin inside the
intron, cutting cox1 in half across position 1. That is the very failure this
rotation exists to prevent, just moved from nad5 to cox1. So the anchor is taken
from the HSP with the lowest qstart (cox1's true 5' exon), and the rotation is
declined altogether when even that HSP starts too far into the protein for the
back-extrapolation to locate the 5' end.

Fail-safe: if no confident cox1 hit is found, or the input is not a single
contig, the assembly is written through UNROTATED (exit 0) so MITOS2 still runs;
mitos_to_emma.py's origin-spanning fallback then handles any residual wrap.

Usage:
    rotate_to_cox1.py --genome asm.fa --cox1-ref cox1_anthozoa.faa --out rotated.fa
"""

import argparse
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# tblastn confidence floor: require a hit covering a decent chunk of cox1 (~505
# aa) at a plausible identity before we trust it enough to re-origin on.
MIN_ALN_AA = 120
MIN_PIDENT = 50.0

# How far into the cox1 protein the anchor HSP may start before back-extrapolating
# to the N-terminus stops being meaningful. A real 5' exon alignment begins within
# the first few residues; anything past this is a 3'-exon-only hit.
MAX_EXTRAPOLATE_AA = 30


def log(msg):
    print(f"[rotate_to_cox1] {msg}", flush=True)


def best_cox1_hit(cox1_ref, genome_fa):
    """tblastn the cox1 protein panel against the genome; pick the anchor HSP.

    Returns (hit, locus, n_exons) where `hit` is a dict
    {sstart, send, qstart, sframe, bitscore, length, pident} in genome (subject)
    nucleotide coordinates, `locus` is the (lo, hi) subject span covered by all
    the anchor query's HSPs, and `n_exons` is how many HSPs that query
    contributed. Returns (None, None, 0) if nothing passes the floor.

    The anchor is the lowest-qstart HSP of the best-scoring *query*, not the
    best-scoring HSP overall: on an intron-split cox1 the 3' exon commonly wins
    on bitscore, and anchoring there puts the origin inside the intron.
    """
    cmd = [
        "tblastn",
        "-query", str(cox1_ref),
        "-subject", str(genome_fa),
        "-seg", "no",
        "-max_target_seqs", "50",
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send sframe bitscore",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log(f"WARNING: tblastn failed (rc={proc.returncode}): {proc.stderr.strip()}")
        return None, None, 0

    by_query = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        hit = {
            "pident": float(f[2]), "length": int(f[3]),
            "qstart": int(f[4]), "qend": int(f[5]),
            "sstart": int(f[6]), "send": int(f[7]),
            "sframe": int(f[8]), "bitscore": float(f[9]),
        }
        if hit["length"] < MIN_ALN_AA or hit["pident"] < MIN_PIDENT:
            continue
        by_query.setdefault(f[0], []).append(hit)
    if not by_query:
        return None, None, 0

    # Pick the query with the most total signal, then keep only its HSPs that
    # agree on strand with its strongest one -- a stray opposite-strand HSP is
    # noise, not an exon.
    qid = max(by_query, key=lambda q: sum(h["bitscore"] for h in by_query[q]))
    hsps = by_query[qid]
    lead = max(hsps, key=lambda h: h["bitscore"])
    plus = lead["sframe"] > 0
    hsps = [h for h in hsps if (h["sframe"] > 0) == plus]

    anchor = min(hsps, key=lambda h: h["qstart"])
    lo = min(min(h["sstart"], h["send"]) for h in hsps)
    hi = max(max(h["sstart"], h["send"]) for h in hsps)
    return anchor, (lo, hi), len(hsps)


def cox1_five_prime(hit, seq_len):
    """Genomic 1-based coordinate of cox1's 5' end and its strand.

    Back-extrapolates from the HSP to the protein N-terminus (qstart) so the
    rotation lands on (or just before) the cox1 start codon rather than mid-gene.
    """
    plus = hit["sframe"] > 0
    extrapolate = 3 * (hit["qstart"] - 1)  # nt upstream of the HSP to the N-term
    if plus:
        pos = hit["sstart"] - extrapolate          # 5' is the lower coord
    else:
        pos = hit["sstart"] + extrapolate          # minus strand: 5' is the higher coord
    # Wrap into 1..seq_len (circular).
    pos = ((pos - 1) % seq_len) + 1
    return pos, plus


def rotate_record(record, start_pos, plus):
    """Return a new SeqRecord re-origined so cox1's 5' end is position 1.

    On a minus-strand cox1 the whole molecule is reverse-complemented first so
    the annotated genome reads with cox1 on the plus strand.
    """
    seq = record.seq
    seq_len = len(seq)
    if not plus:
        seq = seq.reverse_complement()
        start_pos = seq_len - start_pos + 1  # remap coordinate into rc frame
    idx = (start_pos - 1) % seq_len
    rotated = seq[idx:] + seq[:idx]
    return SeqRecord(rotated, id=record.id, description=record.description)


def write_records(records, out_path):
    SeqIO.write(records, str(out_path), "fasta")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, type=Path, help="Assembly FASTA")
    ap.add_argument("--cox1-ref", required=True, type=Path, help="cox1 protein panel (FASTA)")
    ap.add_argument("--out", required=True, type=Path, help="Output (rotated) FASTA")
    args = ap.parse_args()

    if not args.genome.exists():
        sys.exit(f"Genome not found: {args.genome}")
    if not args.cox1_ref.exists():
        sys.exit(f"cox1 reference not found: {args.cox1_ref}")

    records = list(SeqIO.parse(str(args.genome), "fasta"))
    if not records:
        sys.exit(f"No sequences read from {args.genome}")

    # Rotation only makes sense for a single circular molecule.
    if len(records) != 1:
        log(f"{len(records)} records in {args.genome.name}; not a single contig -- "
            "writing through UNROTATED.")
        write_records(records, args.out)
        return

    record = records[0]
    seq_len = len(record.seq)
    hit, locus, n_exons = best_cox1_hit(args.cox1_ref, args.genome)
    if hit is None:
        log("no confident cox1 hit -- writing through UNROTATED (MITOS + "
            "origin-spanning fallback will handle any wrap).")
        write_records([record], args.out)
        return
    if n_exons > 1:
        log(f"cox1 looks intron-split ({n_exons} HSPs, subject {locus[0]}-{locus[1]}); "
            f"anchoring on the 5'-most exon (cox1 residue {hit['qstart']}).")

    # Never cut cox1 in half. start_pos is the anchor HSP's start walked back by
    # 3*(qstart-1) to reach the protein's N-terminus, so it is only trustworthy
    # while that walk is short. A long walk means cox1's real 5' end was never
    # aligned, and the extrapolation is then just as likely to land mid-gene (or
    # in the neighbouring gene) as on the start codon -- which is exactly how an
    # intron-split cox1 used to get re-origined into two pieces.
    if hit["qstart"] - 1 > MAX_EXTRAPOLATE_AA:
        log(f"anchor HSP starts at cox1 residue {hit['qstart']}, too far into the "
            f"protein to locate the 5' end (limit {MAX_EXTRAPOLATE_AA}) -- writing "
            "through UNROTATED rather than risking a rotation inside cox1.")
        write_records([record], args.out)
        return

    start_pos, plus = cox1_five_prime(hit, seq_len)

    log(f"cox1 located: strand={'+' if plus else '-'} 5'pos={start_pos} "
        f"(bitscore={hit['bitscore']:.0f}, {hit['length']}aa, {hit['pident']:.1f}% id); "
        f"re-origining {record.id} ({seq_len} bp).")
    rotated = rotate_record(record, start_pos, plus)
    write_records([rotated], args.out)
    log(f"wrote rotated assembly to {args.out}")


if __name__ == "__main__":
    main()
