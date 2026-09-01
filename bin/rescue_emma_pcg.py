#!/usr/bin/env python3
"""Recover a single EMMA-dropped protein-coding gene (ND4L or ATP8) in place.

EMMA's ``rationalise_matches!`` discards a short CDS when its computed circular
overlap with a longer neighbour exceeds half the shorter feature's length, which
routinely loses ND4L (vs ND4) and ATP8 (vs ATP6). The gene is present in the
assembly; only the annotation is short. This script rebuilds the missing feature
from the flanking-gene coordinates EMMA already produced:

    * define the intergenic window between the two REF-order neighbours,
    * tblastn a small reference-protein set into that window to fix the reading
      frame and get an identity/coverage guard,
    * refine to a clean ORF (vertebrate-mito start + stop codons, with EMMA's
      polyadenylation convention when the stop is completed by the poly-A tail),
    * write matching gene/mRNA/CDS lines into the EMMA .gff and .tbl and the
      per-gene nucleotide/protein FASTAs under cds/ and proteins/.

Every edit is guarded (BLAST identity/coverage, ORF cleanliness, length within
tolerance of the matched reference). A target that fails any guard is left
untouched. A status file records RESCUED / SKIP per target and the script always
exits 0, so a failed rescue reproduces EMMA's original (still-incomplete) bundle
and the assembly continues through LCA / species validation exactly as before.

Runs in the MITOS2 BioContainer (biopython + BLAST+).

Usage:
    rescue_emma_pcg.py --annotation-dir annotation --targets ND4L,ATP8 \\
        --ref-faa rescue_pcg_refs.faa --code 2 \\
        --status <prefix>.rescue.status.txt
"""

import argparse
import re
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from orf_utils import refine_orf

# REF-order neighbours (left, right) whose EMMA coordinates bound the window, and
# the NCBI /product string, for each rescuable gene.
TARGETS = {
    "ND4L": {"left": "TR", "right": "ND4", "product": "NADH dehydrogenase subunit 4L",
             "aa_min": 80, "aa_max": 115, "up_slack": 6, "down_slack": 15},
    "ATP8": {"left": "TK", "right": "ATP6", "product": "ATP synthase F0 subunit 8",
             "aa_min": 40, "aa_max": 80, "up_slack": 6, "down_slack": 48},
}

# Start/stop codon sets and the ORF refiner come from bin/orf_utils.py -- one
# source of truth shared with the annotation QC gate and the coral fixer. This
# script previously carried its own transl_table-2-only tuples while accepting a
# --code argument, so any other table translated under one code and detected
# stops under another; and its own upstream start search kept the FURTHEST
# initiator within 12 codons rather than the nearest, which over-extended the 5'
# end (OG811 ATP8 came out as VKMPQLNP... on a GTG two codons upstream of the
# true ATG). orf_utils.refine_orf ranks candidates properly: a real encoded stop
# first, then canonical ATG over an alternative initiator, then closest to the
# alignment hint.

POLY_A_NOTE = ("putative TAA stop codon is completed by the addition of 3' A "
               "residues to the mRNA")


# --------------------------------------------------------------------------- IO

def find_one(directory, pattern):
    hits = sorted(directory.glob(pattern))
    if not hits:
        raise FileNotFoundError(f"no {pattern} in {directory}")
    return hits[0]


def parse_gff_genes(gff_path):
    """gene name (no MT-) -> (start, end, strand), 1-based inclusive."""
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] != "gene":
                continue
            attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
            name = attrs.get("Name", "").replace("MT-", "")
            if name and name not in genes:
                genes[name] = (int(p[3]), int(p[4]), p[6])
    return genes


def sequence_region(gff_path):
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                return int(line.split()[3])
    return None


# ------------------------------------------------------------------------ BLAST

def tblastn_window(ref_faa, gene, window_fa, code):
    """Best HSP of the gene's reference proteins vs the window nucleotide.

    Returns dict(pid, qcov, ref_id, ref_len, s_from, s_to, frame) or None.
    s_from/s_to are 1-based on the window, 5'->3' along the CDS (s_from < s_to).
    """
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        keep = 0
        for rec in SeqIO.parse(str(ref_faa), "fasta"):
            if rec.id.startswith(f"{gene}_"):
                fh.write(f">{rec.id}\n{rec.seq}\n")
                keep += 1
        qpath = fh.name
    if keep == 0:
        Path(qpath).unlink(missing_ok=True)
        return None

    cols = "qseqid pident length qlen sstart send sframe evalue bitscore"
    proc = subprocess.run(
        ["tblastn", "-query", qpath, "-subject", str(window_fa),
         "-db_gencode", str(code), "-seg", "no", "-max_target_seqs", "50",
         "-evalue", "1e-3", "-outfmt", f"6 {cols}"],
        capture_output=True, text=True)
    Path(qpath).unlink(missing_ok=True)

    best = None
    for line in proc.stdout.strip().splitlines():
        f = line.split("\t")
        if len(f) != 9:
            continue
        pid, aln, qlen = float(f[1]), int(f[2]), int(f[3])
        sstart, send, frame = int(f[4]), int(f[5]), int(f[6])
        evalue, bits = float(f[7]), float(f[8])
        lo, hi = (sstart, send) if sstart <= send else (send, sstart)
        cand = {"pid": pid, "qcov": aln / qlen, "ref_id": f[0], "ref_len": qlen,
                "s_from": lo, "s_to": hi, "frame": frame, "evalue": evalue, "bits": bits}
        if best is None or bits > best["bits"]:
            best = cand
    return best


# -------------------------------------------------------------------- ORF logic
#
# refine_orf lives in bin/orf_utils.py (see the import above): it takes the same
# (window, s_from, down_limit, code) and returns the same
# (cds_start, cds_end, aa, start_codon, poly_a).


# ---------------------------------------------------------------- write patches

def emma_evalue(x):
    return f"{x:.1e}".replace("e-0", "e-").replace("e+0", "e+")


def append_gff(gff_path, chrom, gene, product, g_start, g_end, strand, evalue, poly_a):
    guid, muid, cuid = (uuid.uuid4() for _ in range(3))
    note = f";Note={POLY_A_NOTE}" if poly_a else ""
    with open(gff_path, "a") as out:
        out.write(f"{chrom}\tEmma\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t"
                  f"ID={guid};Name=MT-{gene}\n")
        out.write(f"{chrom}\tEmma\tmRNA\t{g_start}\t{g_end}\t.\t{strand}\t.\t"
                  f"ID={muid};Parent={guid};Name=MT-{gene}\n")
        out.write(f"{chrom}\tEmma\tCDS\t{g_start}\t{g_end}\t{evalue}\t{strand}\t0\t"
                  f"ID={cuid};Parent={muid};Name=MT-{gene};Product={product}{note}\n")


def append_tbl(tbl_path, gene, product, g_start, g_end, strand, code, poly_a):
    a, b = (g_start, g_end) if strand == "+" else (g_end, g_start)
    with open(tbl_path, "a") as out:
        out.write(f"{a}\t{b}\tgene\n\t\t\tgene\tMT-{gene}\n")
        out.write(f"{a}\t{b}\tmRNA\n\t\t\tgene\tMT-{gene}\n")
        out.write(f"{a}\t{b}\tCDS\n")
        out.write(f"\t\t\tproduct\t{product}\n")
        out.write(f"\t\t\ttransl_table\t{code}\n")
        out.write(f"\t\t\tprotein_id\tgnl|Emma|{uuid.uuid4()}\n")
        if poly_a:
            out.write(f"\t\t\tnote\t{POLY_A_NOTE}\n")


def write_fasta(path, header, seq, width=70):
    with open(path, "w") as out:
        out.write(f">{header}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i:i + width] + "\n")


# ----------------------------------------------------------------------- driver

def rescue_gene(gene, genes, genome, chrom, seq_len, ann_dir, ann_name,
                ref_faa, code, min_pid, min_cov, len_tol):
    spec = TARGETS[gene]
    if gene in genes:
        return ("SKIP", gene, "already annotated")
    for flank in (spec["left"], spec["right"]):
        if flank not in genes:
            return ("SKIP", gene, f"flank {flank} absent")

    l_start, l_end, _l_strand = genes[spec["left"]]
    r_start, r_end, r_strand = genes[spec["right"]]
    strand = r_strand  # ND4L/ATP8 sit on the same strand as their downstream PCG

    if strand == "+":
        w_start = l_end + 1 - spec["up_slack"]
        w_end = r_start + spec["down_slack"]
    else:
        w_start = r_end + 1 - spec["down_slack"]
        w_end = l_start + spec["up_slack"]

    if w_start < 1 or w_end > seq_len or w_start >= w_end:
        return ("SKIP", gene, "origin_wrap or degenerate window (re-origin needed)")

    window = genome[w_start - 1:w_end]
    if strand == "-":
        window = str(Seq(window).reverse_complement())
    window = window.upper()

    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(f">window\n{window}\n")
        wpath = fh.name
    try:
        hit = tblastn_window(ref_faa, gene, wpath, code)
    finally:
        Path(wpath).unlink(missing_ok=True)
    if hit is None:
        return ("SKIP", gene, "no tblastn hit in window")
    if hit["pid"] < min_pid or hit["qcov"] < min_cov:
        return ("SKIP", gene,
                f"weak hit pid={hit['pid']:.1f} cov={hit['qcov']:.2f}")

    # Downstream limit in window coordinates: the first base of the downstream
    # gene (plus a codon of the canonical overlap), so ND4L may overlap ND4's
    # start but not run through it.
    if strand == "+":
        down_limit = (r_start - w_start + 1) + 9
    else:
        down_limit = (w_end - r_end + 1) + 9
    down_limit = min(down_limit, len(window))

    orf = refine_orf(window, hit["s_from"], down_limit, code)
    if orf is None:
        return ("SKIP", gene, "no clean ORF in window")
    cds_lo, cds_hi, protein, start_codon, poly_a = orf

    if "*" in protein:
        return ("SKIP", gene, "internal stop in recovered ORF")
    aa_len = len(protein)
    ref_len = hit["ref_len"]
    if not (spec["aa_min"] <= aa_len <= spec["aa_max"]):
        return ("SKIP", gene, f"aa_len {aa_len} outside [{spec['aa_min']},{spec['aa_max']}]")
    if not ((1 - len_tol) * ref_len <= aa_len <= (1 + len_tol) * ref_len):
        return ("SKIP", gene, f"aa_len {aa_len} vs ref {ref_len} (tol {len_tol})")
    if start_codon == "partial":
        return ("SKIP", gene, "5'-partial ORF (no start codon in window)")

    # Map window coords back to genomic 1-based inclusive.
    if strand == "+":
        g_start = w_start + cds_lo - 1
        g_end = w_start + cds_hi - 1
    else:
        g_start = w_end - cds_hi + 1
        g_end = w_end - cds_lo + 1

    cds_nt = genome[g_start - 1:g_end]
    if strand == "-":
        cds_nt = str(Seq(cds_nt).reverse_complement())
    # The gene/CDS feature spans the terminal stop codon, but EMMA's own cds/
    # FASTAs exclude it (every EMMA cds/ file is exactly 3x its protein), so trim
    # it here too -- a rescued gene must not be the one file in the bundle that
    # follows a different convention. A poly-A-completed stop has nothing to trim.
    if not poly_a and len(cds_nt) == 3 * (len(protein) + 1):
        cds_nt = cds_nt[:-3]

    evalue = emma_evalue(hit["evalue"]) if hit["evalue"] > 0 else "0.0"
    append_gff(find_one(ann_dir, "*.gff"), chrom, gene, spec["product"],
               g_start, g_end, strand, evalue, poly_a)
    append_tbl(find_one(ann_dir, "*.tbl"), gene, spec["product"],
               g_start, g_end, strand, code, poly_a)
    (ann_dir / "cds").mkdir(exist_ok=True)
    (ann_dir / "proteins").mkdir(exist_ok=True)
    write_fasta(ann_dir / "cds" / f"MT-{gene}.{ann_name}.fa", ann_name, cds_nt.upper())
    write_fasta(ann_dir / "proteins" / f"MT-{gene}.{ann_name}.fa", ann_name, protein)

    tag = "RESCUED_POLYA" if poly_a else "RESCUED"
    return (tag, gene,
            f"pid={hit['pid']:.1f} cov={hit['qcov']:.2f} ref={hit['ref_id']} "
            f"aa={aa_len} start={start_codon} coords={g_start}..{g_end}({strand})")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True, type=Path)
    ap.add_argument("--targets", required=True, help="comma list from {ND4L,ATP8}")
    ap.add_argument("--ref-faa", required=True, type=Path)
    ap.add_argument("--code", type=int, default=2)
    ap.add_argument("--status", required=True, type=Path)
    ap.add_argument("--min-pid", type=float, default=55.0)
    ap.add_argument("--min-cov", type=float, default=0.75)
    ap.add_argument("--len-tol", type=float, default=0.15)
    args = ap.parse_args()

    ann_dir = args.annotation_dir
    lines = []

    def finish():
        args.status.parent.mkdir(parents=True, exist_ok=True)
        args.status.write_text("".join(f"{s}\t{g}\t{m}\n" for s, g, m in lines))
        for s, g, m in lines:
            print(f"[rescue_emma_pcg] {s}\t{g}\t{m}", file=sys.stderr)
        sys.exit(0)

    try:
        gff = find_one(ann_dir, "*.gff")
    except FileNotFoundError as exc:
        lines.append(("SKIP", "-", f"no GFF: {exc}"))
        finish()

    ann_name = gff.stem  # e.g. OG624.ilmn.230419.getorg1770.emma102
    try:
        genes = parse_gff_genes(gff)
        with open(find_one(ann_dir, "*.fa")) as fh:
            rec = next(SeqIO.parse(fh, "fasta"))
            genome = str(rec.seq).upper()
            chrom = rec.id
        seq_len = sequence_region(gff) or len(genome)
    except Exception as exc:  # noqa: BLE001 - degrade to no-op
        lines.append(("SKIP", "-", f"parse error: {exc}"))
        finish()

    for gene in [t.strip().upper() for t in args.targets.split(",") if t.strip()]:
        if gene not in TARGETS:
            lines.append(("SKIP", gene, "not a rescuable gene"))
            continue
        try:
            lines.append(rescue_gene(gene, genes, genome, chrom, seq_len, ann_dir,
                                     ann_name, args.ref_faa, args.code,
                                     args.min_pid, args.min_cov, args.len_tol))
        except Exception as exc:  # noqa: BLE001 - one target's failure is a SKIP
            lines.append(("SKIP", gene, f"error: {exc}"))

    if not lines:
        lines.append(("SKIP", "-", "no targets"))
    finish()


if __name__ == "__main__":
    main()
