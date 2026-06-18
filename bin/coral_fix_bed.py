#!/usr/bin/env python3
"""
Repair an anthozoan (coral) MITOS2 result.bed before it is handed to
mitos_to_emma.py.

MITOS2 reliably annotates the 13 PCGs + 12S of coral mitogenomes but, on
divergent anthozoans, routinely (a) fails to call the 16S rRNA (rrnL) and
(b) reports only one exon of the group-I-intron-split nad5. Both features are,
however, present in the assembly and well conserved in close coral references.

This script transfers just those two features from a curated close coral
reference (GenBank) by BLAST and rewrites result.bed:

  * inject an ``rrnL`` row (mitos_to_emma maps rrnL -> RNR2 = 16S)
  * replace the partial nad5 row(s) with one row per exon, named nad5_0,
    nad5_1, ... in transcript order, so mitos_to_emma concatenates them into
    the correct spliced CDS + translation
  * drop a spurious tRNA wholly inside the new 16S span (MITOS sometimes mis-
    calls a weak trnV there)

Everything stays in the BED coordinate frame MITOS used (the cox1-rotated
genome). mitos_to_emma.py then does the gff / cds/ / proteins/ standardisation
and the final re-origin to trnM, unchanged.

Fail-safe: QC gates (BLAST coverage/identity; reconstructed nad5 must be a clean
ORF under the genetic code) guard every edit. If a gate fails the corresponding
feature is left untouched; the original rows are never corrupted. A status line
(FIXED / PARTIAL / SKIPPED / FAILED + reason) is written to --status, and the
script always exits 0 so a poor reference can never break the run.

Runs in the MITOS2 BioContainer (provides blastn + biopython).

Usage:
    coral_fix_bed.py --bed result.bed --genome cox1_rotated.fa \
        --ref-gb ref.reference.gb --code 4 \
        --out-bed result.fixed.bed --status coral_fix.status.txt
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


# ---- reference parsing -----------------------------------------------------

def ref_features(gb_path):
    """Return (rrnl_seq, nad5_exon_seqs) extracted from the reference GenBank.

    nad5_exon_seqs is a list of exon nucleotide strings in transcript order.
    Either may be None/empty if the reference lacks that (well-)annotated feature.
    """
    rec = next(SeqIO.parse(str(gb_path), "genbank"))
    full = str(rec.seq).upper()

    def label(feat):
        q = feat.qualifiers
        return " ".join(q.get("gene", []) + q.get("product", [])).upper()

    rrnl = None
    nad5_exons = None
    for feat in rec.features:
        lab = label(feat)
        if feat.type == "rRNA" and ("16S" in lab or "RRNL" in lab or "LARGE" in lab):
            rrnl = str(feat.extract(rec.seq)).upper()
        if feat.type == "CDS" and ("ND5" in lab or "NAD5" in lab or "SUBUNIT 5" in lab):
            # one entry per exon (a coral nad5 has 2), already in transcript order
            exons = []
            for part in feat.location.parts:
                seg = full[int(part.start):int(part.end)]
                if part.strand == -1:
                    seg = str(Seq(seg).reverse_complement())
                exons.append(seg)
            nad5_exons = exons
    return rrnl, nad5_exons


# ---- BLAST -----------------------------------------------------------------

def blast_probes(probes, genome_fa):
    """probes: dict name->seq. Returns name->(pid, sstart, send, scov) best hit."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for k, v in probes.items():
            fh.write(f">{k}\n{v}\n")
        qpath = fh.name
    out = subprocess.run(
        ["blastn", "-query", qpath, "-subject", str(genome_fa),
         "-evalue", "1e-5", "-word_size", "9",
         "-outfmt", "6 qseqid pident length sstart send qlen"],
        capture_output=True, text=True).stdout
    Path(qpath).unlink(missing_ok=True)
    best = {}
    for line in out.strip().splitlines():
        q, pid, ln, ss, se, qlen = line.split("\t")
        ln = int(ln)
        if q not in best or ln > best[q][0]:
            best[q] = (ln, float(pid), int(ss), int(se), ln / int(qlen))
    # drop the sort length, return (pid, sstart, send, scov)
    return {k: (v[1], v[2], v[3], v[4]) for k, v in best.items()}


# ---- translation -----------------------------------------------------------

def clean_orf(genome_seq, exons_1based, code):
    """exons_1based: list of (start,end,strand) 1-based inclusive in transcript
    order. Reconstruct, translate, and report (ok, length_nt, aa, n_internal_stops)."""
    nt = ""
    for s, e, st in exons_1based:
        seg = genome_seq[s - 1:e]
        if st == "-":
            seg = str(Seq(seg).reverse_complement())
        nt += seg
    aa = str(Seq(nt).translate(table=code))
    internal = aa[:-1].count("*")
    ends_stop = aa.endswith("*")
    return (internal == 0 and ends_stop), len(nt), aa, internal


# ---- BED helpers -----------------------------------------------------------

TRNA_RE = re.compile(r"^trn", re.I)

def parse_bed(path):
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cols = line.split("\t")
        if len(cols) < 6:
            cols = line.split()
        if len(cols) < 6:
            continue
        rows.append(cols[:6])  # chrom,start0,end,name,score,strand
    return rows

def bare_name(name):
    return re.sub(r"[-_](?:\d+|[a-z])$", "", re.sub(r"\(.*?\)$", "", name.strip()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--genome", required=True, type=Path, help="cox1-rotated MITOS input fasta")
    ap.add_argument("--ref-gb", required=True, type=Path)
    ap.add_argument("--code", type=int, default=4)
    ap.add_argument("--out-bed", required=True, type=Path)
    ap.add_argument("--status", required=True, type=Path)
    ap.add_argument("--min-cov", type=float, default=0.80)
    ap.add_argument("--min-pid", type=float, default=80.0)
    ap.add_argument("--nad5-min-gain", type=int, default=30,
                    help="Replace MITOS's present nad5 with the reference-transferred "
                         "join when the rebuilt clean ORF is at least this many nt "
                         "longer (catches a truncated but present first exon).")
    args = ap.parse_args()

    def finish(state, msg, rows):
        with open(args.out_bed, "w") as out:
            for r in rows:
                out.write("\t".join(map(str, r)) + "\n")
        args.status.write_text(f"{state}\t{msg}\n")
        print(f"[coral_fix] {state}: {msg}", file=sys.stderr)
        sys.exit(0)

    rows = parse_bed(args.bed)
    if not rows:
        finish("FAILED", "empty or unreadable result.bed", rows)
    chrom = rows[0][0]

    genome_rec = next(SeqIO.parse(str(args.genome), "fasta"))
    gseq = str(genome_rec.seq).upper()

    try:
        rrnl_seq, nad5_exons = ref_features(args.ref_gb)
    except Exception as exc:  # malformed reference -> leave MITOS output as-is
        finish("FAILED", f"could not parse reference {args.ref_gb.name}: {exc}", rows)

    probes = {}
    if rrnl_seq:
        probes["rrnL"] = rrnl_seq
    if nad5_exons:
        for i, ex in enumerate(nad5_exons):
            probes[f"nad5_{i}"] = ex
    if not probes:
        finish("FAILED", "reference lacks usable 16S and nad5 annotation", rows)

    hits = blast_probes(probes, args.genome)

    have_rrnl = any(bare_name(r[3]) == "rrnL" for r in rows)
    nad5_rows = [r for r in rows if bare_name(r[3]) == "nad5"]

    actions = []
    new_rrnl_span = None

    # ---- 16S ----
    add_rrnl = None
    if not have_rrnl and "rrnL" in hits:
        pid, ss, se, cov = hits["rrnL"]
        if cov >= args.min_cov and pid >= args.min_pid:
            lo, hi = min(ss, se), max(ss, se)
            strand = "+" if ss < se else "-"
            add_rrnl = [chrom, lo - 1, hi, "rrnL", f"{pid:.1f}", strand]
            new_rrnl_span = (lo, hi)
            actions.append(f"+16S {lo}-{hi}({strand}) cov={cov:.2f} pid={pid:.1f}")
        else:
            actions.append(f"16S-skip cov={cov:.2f} pid={pid:.1f}")

    # ---- nad5 ----
    add_nad5 = None
    if nad5_exons:
        exon_hits = [hits.get(f"nad5_{i}") for i in range(len(nad5_exons))]
        if all(h is not None for h in exon_hits):
            spans = []
            ok_thresh = True
            for h in exon_hits:
                pid, ss, se, cov = h
                if cov < args.min_cov or pid < args.min_pid:
                    ok_thresh = False
                spans.append((min(ss, se), max(ss, se), "+" if ss < se else "-"))
            if ok_thresh:
                exons_1b = [(s, e, st) for (s, e, st) in spans]  # transcript order = ref order
                ok, ln, aa, internal = clean_orf(gseq, exons_1b, args.code)
                # Rebuild when MITOS is missing an exon (fewer rows than the
                # reference) OR when MITOS has all exons but a truncated one: the
                # group-I-intron first exon is routinely under-called, leaving a
                # present-but-short nad5 the gate can still wave through. Compare
                # the rebuilt clean ORF against MITOS's current nad5 nt length and
                # take the reference join when it is meaningfully longer.
                cur_nt = sum(int(r[2]) - int(r[1]) for r in nad5_rows)  # bed: end - start0
                missing_exon = len(nad5_rows) < len(nad5_exons)
                longer = ok and cur_nt and (ln - cur_nt) >= args.nad5_min_gain
                if ok and (missing_exon or longer or not nad5_rows):
                    add_nad5 = []
                    for i, (s, e, st) in enumerate(spans):
                        add_nad5.append([chrom, s - 1, e, f"nad5_{i}", "0.0", st])
                    why = ("missing-exon" if missing_exon
                           else f"longer-orf {cur_nt}->{ln}nt" if longer else "rebuilt")
                    actions.append(f"nad5 join {'+'.join(f'{s}-{e}' for s,e,_ in spans)} "
                                   f"{len(aa)}aa stops={internal} ({why})")
                elif missing_exon and not ok:
                    actions.append(f"nad5-skip internal_stops={internal} (kept MITOS nad5)")
                else:
                    actions.append(f"nad5-ok (MITOS {cur_nt}nt; "
                                   f"rebuilt {ln if ok else 'NA'}nt clean={ok})")
            else:
                actions.append("nad5-skip below thresholds")

    if add_rrnl is None and add_nad5 is None:
        finish("SKIPPED", ("; ".join(actions)) or "nothing to add", rows)

    # ---- build patched bed ----
    out_rows = []
    for r in rows:
        bn = bare_name(r[3])
        if add_nad5 is not None and bn == "nad5":
            continue  # replaced below
        if add_rrnl is not None and new_rrnl_span and TRNA_RE.match(r[3]):
            s0, e1 = int(r[1]) + 1, int(r[2])
            if s0 >= new_rrnl_span[0] and e1 <= new_rrnl_span[1]:
                actions.append(f"-{r[3]} (spurious tRNA inside 16S)")
                continue
        out_rows.append(r)
    if add_rrnl is not None:
        out_rows.append(add_rrnl)
    if add_nad5 is not None:
        out_rows.extend(add_nad5)

    state = "FIXED" if (add_rrnl is not None and add_nad5 is not None) else "PARTIAL"
    finish(state, "; ".join(actions), out_rows)


if __name__ == "__main__":
    main()
