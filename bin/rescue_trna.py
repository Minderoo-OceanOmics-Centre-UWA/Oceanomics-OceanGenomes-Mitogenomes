#!/usr/bin/env python3
"""Recover EMMA-dropped tRNAs in place with a second, independent tRNA finder.

EMMA's covariance model periodically misses a tRNA that is physically present in
the assembly on an otherwise complete, correctly ordered vertebrate mitogenome.
The gate (bin/trna_rescue_gate.py) only routes an assembly here when every missing
REF gene is a tRNA and the whole 13-PCG + 2-rRNA core is present and ordered, so
the insertion gap for each target is well defined.

tRNAscan-SE 2.0 (vertebrate-mitochondrial model) is run against the EMMA genome
FASTA by the upstream TRNA_SCAN process -- its Perl BioContainer has no Python, so
the scan and this stdlib-only splicer run in separate containers. This script
reads that scan's tabular output (--scan-out) and, for each requested missing
tRNA, splices back the single best hit that survives every guard:

    * isotype AND anticodon match the specific missing gene (disambiguates the
      two Leu and two Ser isotypes),
    * Infernal score >= --min-score,
    * length in [--min-len, --max-len] and intronless,
    * the hit's midpoint lies in the genomic gap between the target's nearest
      present REF-order neighbours (origin-spanning gap -> skip),
    * the hit overlaps no existing annotated feature by more than --max-overlap bp,
    * exactly one unambiguous surviving hit.

A surviving hit is written as gene + tRNA lines into the EMMA .gff and a gene +
tRNA block into the .tbl, matching how EMMA writes its own tRNAs. No cds/ or
proteins/ FASTA is produced (tRNAs have none). Every target that fails a guard is
left untouched; a status file records RESCUED / SKIP per target and the script
always exits 0, so a failed rescue reproduces EMMA's original (still-incomplete)
bundle and the assembly is held at the QC gate exactly as before -- unless the
residual shortfall is now within annotation_trna_tolerance, in which case
annotation_stats.py lets it pass with the residue recorded in trna_advisory.

Stdlib only -- runs in the psycopg2 BioContainer alongside annotation_stats.py.

Usage:
    rescue_trna.py --annotation-dir annotation --targets TW,TA,TN \\
        --scan-out <prefix>.trnascan.tsv --status <prefix>.trna_rescue.status.txt
"""

import argparse
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path

# REF tRNA name -> (tRNAscan-SE isotype 'Type', accepted anticodons in DNA), and
# the canonical /product string, both derived from bin/mito_gene_order.py so the
# products written here are byte-identical to the rest of the pipeline's.
# REF_GENES is shared with the gates and annotation_stats.py for the same reason.
from mito_gene_order import (
    REF_GENES,
    TRNA_AA,
    TRNA_ANTICODON,
    TRNA_PRODUCT as PRODUCT,
)

TRNA_SPEC = {
    name: (TRNA_AA[name[1:]], {ac}) for name, ac in TRNA_ANTICODON.items()
}


FEATURE_TYPES = {"gene", "cds", "trna", "rrna"}


# --------------------------------------------------------------------------- IO

def find_one(directory, pattern):
    hits = sorted(directory.glob(pattern))
    if not hits:
        raise FileNotFoundError(f"no {pattern} in {directory}")
    return hits[0]


def read_single_fasta(path):
    """Return (header_id, sequence) for a single-record FASTA. stdlib only."""
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].split()[0]
            elif header is not None:
                seq.append(line.strip())
    return header, "".join(seq).upper()


def parse_gff(gff_path):
    """(genes, features) where genes maps REF gene name -> (start,end,strand)
    1-based inclusive, and features is a list of (lo,hi) intervals for every
    gene/CDS/tRNA/rRNA feature line."""
    genes = {}
    features = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2].lower() not in FEATURE_TYPES:
                continue
            start, end = int(p[3]), int(p[4])
            features.append((min(start, end), max(start, end)))
            if p[2] != "gene":
                continue
            attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
            name = attrs.get("Name", "").replace("MT-", "")
            if name and name not in genes:
                genes[name] = (start, end, p[6])
    return genes, features


# ------------------------------------------------------------------- tRNAscan-SE

def parse_trnascan(text):
    """Parse tRNAscan-SE tabular (-o) output into a list of hit dicts:
    dict(type, anticodon, lo, hi, strand, score, intron, note); lo/hi are 1-based
    inclusive on the plus strand."""
    hits = []
    for row in text.splitlines():
        f = row.split()
        # data rows start with the sequence name then an integer tRNA number;
        # the 3 header lines and the dashed separator do not.
        if len(f) < 9 or not f[1].isdigit():
            continue
        try:
            begin, end = int(f[2]), int(f[3])
            score = float(f[8])
            intron = int(f[6]) or int(f[7])
        except ValueError:
            continue
        note = " ".join(f[9:]).lower() if len(f) > 9 else ""
        hits.append({
            "type": f[4],
            "anticodon": f[5].upper(),
            "lo": min(begin, end),
            "hi": max(begin, end),
            "strand": "+" if begin <= end else "-",
            "score": score,
            "intron": intron,
            "note": note,
        })
    return hits


def run_trnascan(genome_fa, model):
    """Run tRNAscan-SE with the given mito model and return parsed hit dicts."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "trnascan.out"
        cmd = ["tRNAscan-SE", "-M", model, "-q", "-Q",
               "-o", str(out), str(genome_fa)]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return parse_trnascan(out.read_text())


# ---------------------------------------------------------------------- placement

def neighbour_gap(target, genes):
    """Genomic (lo, hi) between target's nearest present REF-order neighbours,
    or None if a neighbour is missing on either side or the gap wraps the origin."""
    i = REF_GENES.index(target)
    left = next((g for g in reversed(REF_GENES[:i]) if g in genes), None)
    right = next((g for g in REF_GENES[i + 1:] if g in genes), None)
    if left is None or right is None:
        return None
    left_hi = max(genes[left][0], genes[left][1])
    right_lo = min(genes[right][0], genes[right][1])
    if not left_hi < right_lo:
        return None  # origin-spanning window
    return left_hi, right_lo


def overlap_bp(a_lo, a_hi, b_lo, b_hi):
    return max(0, min(a_hi, b_hi) - max(a_lo, b_lo) + 1)


# ---------------------------------------------------------------- write patches

def append_gff(gff_path, chrom, name, product, start, end, strand):
    guid, tuid = uuid.uuid4(), uuid.uuid4()
    with open(gff_path, "a") as out:
        out.write(f"{chrom}\tEmma\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                  f"ID={guid};Name=MT-{name}\n")
        out.write(f"{chrom}\tEmma\ttRNA\t{start}\t{end}\t.\t{strand}\t.\t"
                  f"ID={tuid};Parent={guid};Name=MT-{name};Product={product}\n")


def append_tbl(tbl_path, name, product, start, end, strand, score):
    a, b = (start, end) if strand == "+" else (end, start)
    with open(tbl_path, "a") as out:
        out.write(f"{a}\t{b}\tgene\n\t\t\tgene\tMT-{name}\n")
        out.write(f"{a}\t{b}\ttRNA\n\t\t\tproduct\t{product}\n")
        out.write(f"\t\t\tnote\trecovered by tRNAscan-SE (Infernal score {score:.1f})\n")


# ----------------------------------------------------------------------- driver

def rescue_one(target, hits, genes, features, args):
    if target not in TRNA_SPEC:
        return ("SKIP", target, "not a known mt tRNA")
    if target in genes:
        return ("SKIP", target, "already annotated")

    aa, anticodons = TRNA_SPEC[target]
    gap = neighbour_gap(target, genes)
    if gap is None:
        return ("SKIP", target, "neighbour missing or origin-spanning gap")
    gap_lo, gap_hi = gap

    cands = []
    for h in hits:
        if h["type"] != aa or h["anticodon"] not in anticodons:
            continue
        if "pseudo" in h["note"] or "trunc" in h["note"]:
            continue
        if h["intron"]:
            continue
        length = h["hi"] - h["lo"] + 1
        if not (args.min_len <= length <= args.max_len):
            continue
        if h["score"] < args.min_score:
            continue
        mid = (h["lo"] + h["hi"]) / 2
        if not (gap_lo <= mid <= gap_hi):
            continue
        clash = max((overlap_bp(h["lo"], h["hi"], f_lo, f_hi)
                     for f_lo, f_hi in features), default=0)
        if clash > args.max_overlap:
            continue
        cands.append(h)

    if not cands:
        return ("SKIP", target, "no tRNAscan hit survived the guards")

    cands.sort(key=lambda h: h["score"], reverse=True)
    if len(cands) > 1 and (cands[0]["score"] - cands[1]["score"]) < 2.0 \
            and overlap_bp(cands[0]["lo"], cands[0]["hi"],
                           cands[1]["lo"], cands[1]["hi"]) == 0:
        return ("SKIP", target, "ambiguous location (two comparable hits)")

    h = cands[0]
    chrom = args.chrom
    append_gff(find_one(args.ann_dir, "*.gff"), chrom, target, PRODUCT[target],
               h["lo"], h["hi"], h["strand"])
    append_tbl(find_one(args.ann_dir, "*.tbl"), target, PRODUCT[target],
               h["lo"], h["hi"], h["strand"], h["score"])
    # Record the placement in the caller's own state before the next target is
    # considered: without this a second target's neighbour_gap still sees the
    # pre-rescue gene set, and its --max-overlap guard cannot see the tRNA that
    # was just written, so two rescued tRNAs could be placed on top of each other.
    genes[target] = (h["lo"], h["hi"], h["strand"])
    features.append((h["lo"], h["hi"]))
    return ("RESCUED", target,
            f"score={h['score']:.1f} anticodon={h['anticodon']} "
            f"coords={h['lo']}..{h['hi']}({h['strand']})")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True, type=Path)
    ap.add_argument("--targets", required=True, help="comma list of REF tRNA names")
    ap.add_argument("--status", required=True, type=Path)
    ap.add_argument("--scan-out", type=Path,
                    help="tRNAscan-SE tabular (-o) output from the TRNA_SCAN "
                         "process. If omitted, tRNAscan-SE is run directly "
                         "(needs it on PATH).")
    ap.add_argument("--model", default="vert", help="tRNAscan-SE -M model")
    ap.add_argument("--min-score", type=float, default=20.0)
    ap.add_argument("--max-overlap", type=int, default=0)
    ap.add_argument("--min-len", type=int, default=30)
    ap.add_argument("--max-len", type=int, default=100)
    ap.add_argument("--code", type=int, default=2,
                    help="accepted for call-site symmetry; -M vert is code-independent. "
                         "Only the vertebrate (code 2) EMMA path reaches this rescue, so "
                         "2 is the path's only value, not a guess.")
    args = ap.parse_args()

    args.ann_dir = args.annotation_dir
    lines = []

    def finish():
        args.status.parent.mkdir(parents=True, exist_ok=True)
        args.status.write_text("".join(f"{s}\t{g}\t{m}\n" for s, g, m in lines))
        for s, g, m in lines:
            print(f"[rescue_trna] {s}\t{g}\t{m}", file=sys.stderr)
        sys.exit(0)

    try:
        gff = find_one(args.ann_dir, "*.gff")
        genome_fa = find_one(args.ann_dir, "*.fa")
        genes, features = parse_gff(gff)
        args.chrom, _genome = read_single_fasta(genome_fa)
    except Exception as exc:  # noqa: BLE001 - degrade to no-op
        lines.append(("SKIP", "-", f"parse error: {exc}"))
        finish()

    targets = [t.strip().upper() for t in args.targets.split(",") if t.strip()]
    try:
        if args.scan_out:
            hits = parse_trnascan(args.scan_out.read_text())
        else:
            hits = run_trnascan(genome_fa, args.model)
    except Exception as exc:  # noqa: BLE001 - a missing/empty scan is a SKIP for every target
        for t in targets or ["-"]:
            lines.append(("SKIP", t, f"tRNAscan-SE scan unavailable: {exc}"))
        finish()

    for target in targets:
        try:
            lines.append(rescue_one(target, hits, genes, features, args))
        except Exception as exc:  # noqa: BLE001 - one target's failure is a SKIP
            lines.append(("SKIP", target, f"error: {exc}"))

    if not lines:
        lines.append(("SKIP", "-", "no targets"))
    finish()


if __name__ == "__main__":
    main()
