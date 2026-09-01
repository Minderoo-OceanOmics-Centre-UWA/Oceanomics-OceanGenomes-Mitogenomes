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
  * rebuild an intron-split cox1 as cox1_0/cox1_1. Some scleractinians carry a
    second group-I intron in cox1, holding a LAGLIDADG homing endonuclease ORF
    (MITOS calls that ORF `lagli`, correctly). MITOS then reports only cox1's 3'
    exon under the cox1 name -- roughly half the gene -- and misses the 5' exon
    on the far side of the intron, so the CO1 barcode and MT-CO1 protein come
    out truncated. The reference's single-exon cox1 BLASTs onto such an assembly
    as two subject blocks split by the intron, which is what the join needs.

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

from orf_utils import (classify_cds, refine_orf, revcomp, start_codons,
                       stop_codons)

# EMMA-contract PCG bare names in the MITOS BED (nad5 is owned by the
# reference-transfer pass above, so it is excluded from the ORF-snap pass).
PCG_BARE = {
    "cox1", "cox2", "cox3", "cob", "atp6", "atp8",
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad6",
}


# ---- reference parsing -----------------------------------------------------

def ref_features(gb_path):
    """Return (rrnl_seq, nad5_exon_seqs, cox1_seq) from the reference GenBank.

    nad5_exon_seqs is a list of exon nucleotide strings in transcript order.
    cox1_seq is the reference cox1 CDS as one string -- references are normally
    single-exon there, and that is precisely what makes it a usable probe for an
    intron-split cox1 in the sample (see cox1_intron_join). Any of the three may
    be None/empty if the reference lacks that (well-)annotated feature.
    """
    rec = next(SeqIO.parse(str(gb_path), "genbank"))
    full = str(rec.seq).upper()

    def label(feat):
        q = feat.qualifiers
        return " ".join(q.get("gene", []) + q.get("product", [])).upper()

    rrnl = None
    nad5_exons = None
    cox1 = None
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
        # "SUBUNIT I" must not also match subunits II/III.
        if feat.type == "CDS" and ("COX1" in lab or "CO1" in lab
                                   or re.search(r"SUBUNIT I(?![IV])", lab)):
            cox1 = str(feat.extract(rec.seq)).upper()
    return rrnl, nad5_exons, cox1


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


def blast_hsps(name, probe, genome_fa, min_len=60):
    """All HSPs for one probe, not just the best, in query order.

    blast_probes() keeps one hit per query, which is the right call for a probe
    expected to land once. An intron-split gene lands as several blocks and the
    whole point is to see all of them, so this returns every HSP at or above
    `min_len` aligned nt as dicts {qs, qe, ss, se, pid, ln}.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(f">{name}\n{probe}\n")
        qpath = fh.name
    out = subprocess.run(
        ["blastn", "-query", qpath, "-subject", str(genome_fa),
         "-evalue", "1e-5", "-word_size", "9",
         "-outfmt", "6 qstart qend sstart send pident length"],
        capture_output=True, text=True).stdout
    Path(qpath).unlink(missing_ok=True)
    hsps = []
    for line in out.strip().splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        ln = int(f[5])
        if ln < min_len:
            continue
        hsps.append({"qs": int(f[0]), "qe": int(f[1]), "ss": int(f[2]),
                     "se": int(f[3]), "pid": float(f[4]), "ln": ln})
    return hsps


def group_exons(hsps, n, slop=15):
    """Collapse plus-strand HSPs into exon blocks in doubled-genome coordinates.

    Two HSPs belong to the same exon when they sit on the same diagonal: the
    subject gap between them matches the query gap. A group-I intron breaks that
    diagonal (large subject gap, no query gap) and so starts a new exon. The
    subject is also tried one genome-length on, so an exon that runs off the end
    of the contig and continues past position 1 stays a single block instead of
    being mistaken for two exons -- which is the normal case here, because the
    linearisation point routinely falls inside cox1.

    Returns a list of (start, end) 1-based inclusive spans in transcript order,
    possibly exceeding n where a span wraps the origin.
    """
    plus = sorted((h for h in hsps if h["se"] > h["ss"]), key=lambda h: h["qs"])
    kept = []
    for h in plus:
        # Drop an HSP whose query range is already covered: a short spurious hit
        # elsewhere in the genome must not open a third "exon".
        if any(h["qs"] >= k["qs"] and h["qe"] <= k["qe"] for k in kept):
            continue
        kept.append(h)
    exons, prev_qe = [], None
    for h in kept:
        if not exons:
            exons.append([h["ss"], h["se"]])
            prev_qe = h["qe"]
            continue
        qgap = h["qs"] - prev_qe - 1
        for cand in (h["ss"], h["ss"] + n):
            if abs((cand - exons[-1][1] - 1) - qgap) <= slop:
                exons[-1][1] = cand + (h["se"] - h["ss"])
                break
        else:
            exons.append([h["ss"], h["se"]])
        prev_qe = h["qe"]
    return [(s, e) for s, e in exons]


def cox1_intron_join(cox1_ref, genome_fa, g2, n, code,
                     min_cov=0.80, min_pid=80.0, max_scan_codons=60):
    """Rebuild a group-I-intron-split cox1 from the reference cox1 probe.

    Returns (exons, aa, note) with `exons` as 1-based inclusive (start, end)
    spans in doubled-genome coordinates, or (None, None, reason).

    The BLAST blocks give the splice junctions but not the CDS ends: the probe's
    own start codon is transferred with the alignment, while its 3' end stops at
    the last aligned base, which is short of the sample's stop codon. So the 5'
    end is only slid if it did not land on a start codon, and the 3' end is
    walked forward in frame to the first stop.
    """
    hsps = blast_hsps("cox1_ref", cox1_ref, genome_fa)
    if not hsps:
        return None, None, "no reference hit"
    exons = group_exons(hsps, n)
    if len(exons) < 2:
        return None, None, f"not intron-split ({len(exons)} block)"
    if len(exons) > 2:
        return None, None, f"{len(exons)} blocks (expected 2)"

    aligned = sum(h["ln"] for h in hsps if h["se"] > h["ss"])
    cov = min(aligned, len(cox1_ref)) / len(cox1_ref)
    pid = max(h["pid"] for h in hsps)
    if cov < min_cov or pid < min_pid:
        return None, None, f"below thresholds cov={cov:.2f} pid={pid:.1f}"

    def spliced(exs):
        return "".join(g2[s - 1:e] for s, e in exs)

    exons = [list(e) for e in exons]
    # 5': accept the transferred start when it already is one, else slide.
    if g2[exons[0][0] - 1:exons[0][0] + 2] not in start_codons(code):
        for d in (3, -3, 6, -6, 9, -9, 12, -12):
            cand = exons[0][0] + d
            if cand > 0 and g2[cand - 1:cand + 2] in start_codons(code):
                exons[0][0] = cand
                break
        else:
            return None, None, "no start codon near the transferred 5' end"
    # 3': restore frame, then walk to the first in-frame stop.
    end = exons[-1][1] + (3 - sum(e - s + 1 for s, e in exons) % 3) % 3
    for _ in range(max_scan_codons):
        if g2[end - 3:end] in stop_codons(code):
            break
        end += 3
    else:
        return None, None, "no in-frame stop within scan window"
    exons[-1][1] = end

    nt = spliced(exons)
    info = classify_cds(nt, code)
    if not (info["start_ok"] and info["stop_ok"]) or info["internal_stops"]:
        return None, None, (f"rebuilt ORF not clean "
                            f"(internal_stops={info['internal_stops']})")
    return [tuple(e) for e in exons], nt, f"cov={cov:.2f} pid={pid:.1f}"


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


def cds_from_row(gseq, start0, end, strand):
    """Spliced CDS nucleotides for one BED row, handling the circular origin.

    MITOS writes an origin-spanning feature with start > end (e.g. `17684 19` on
    a 17.7 kb coral): the CDS is the tail of the sequence followed by the head.
    """
    n = len(gseq)
    seg = gseq[start0:end] if start0 < end else gseq[start0:] + gseq[:end]
    return revcomp(seg) if strand == "-" else seg.upper()


def _linear_span(start0, end, n):
    """A BED (start0, end) mapped into a doubled-genome coordinate space, where
    an origin-spanning feature becomes a plain interval (start0, end + n)."""
    return (start0, end) if start0 < end else (start0, end + n)


# A PCG that shares a canonical overlap with its neighbour (ATP8 into ATP6,
# ND4L into ND4) must be allowed to snap into that neighbour rather than being
# clipped at its first base, so the neighbour bound is relaxed by this many nt
# for those pairs. Anything more would let a boundary error walk into the
# neighbouring gene.
CANONICAL_OVERLAP = {"atp8": 12, "nad4l": 12}

# Per-gene absolute floor (nt) on how far a snap may move the CDS length. A flat
# floor of 60 nt is over a third of atp8 (~200 nt) and a fifth of nad3/nad4l, so
# short PCGs get a proportionate one instead.
LEN_TOL_FLOOR = {"atp8": 18, "nad4l": 18, "nad3": 24}
DEFAULT_LEN_TOL_FLOOR = 60


def snap_pcg_rows(rows, gseq, code, window=36, len_tol=0.15, skip=()):
    """Per-PCG ORF-boundary snap, assembly-only (no reference).

    For every single-exon PCG row whose reconstructed CDS is not a clean ORF
    under `code`, build a small window bounded by the neighbouring rows and call
    orf_utils.refine_orf. When it yields a clean ORF within `len_tol` of the
    original length, return the rewritten (start0, end) for that row.

    All coordinate work happens in a doubled-genome space so a feature spanning
    the circular origin (MITOS writes those with start > end) is an ordinary
    interval; results are mapped back, restoring the start > end form if the
    snapped CDS still wraps. Before this, every origin-spanning PCG was reported
    `no window` and silently went unchecked -- which in practice meant atp8 was
    never checked on any coral, because ROTATE_ORIGIN puts cox1 at position 1 and
    leaves atp8 straddling the join.

    Bare names in `skip` are left alone: a gene rebuilt by a reference-transfer
    join owns its own coordinates, and snapping the stale single-exon row it
    replaces would only produce a contradictory second opinion.

    Returns (fixes, actions, unrepaired):
        fixes      = {row_index: (new_start0, new_end)}
        actions    = human-readable strings for the status file
        unrepaired = bare names of PCGs found broken that could NOT be repaired,
                     so the caller does not report the sample as FIXED
    """
    fixes = {}
    actions = []
    unrepaired = []

    n = len(gseq)
    gseq2 = (gseq + gseq).upper()

    # Row spans in doubled coordinates, for neighbour clipping. Each row is also
    # considered one genome-length later, so a neighbour on the far side of the
    # origin still bounds a wrapped feature.
    spans = []
    for r in rows:
        try:
            lo, hi = _linear_span(int(r[1]), int(r[2]), n)
        except (ValueError, IndexError):
            continue
        spans.append((lo, hi))
        spans.append((lo + n, hi + n))

    by_bare = {}
    for i, r in enumerate(rows):
        by_bare.setdefault(bare_name(r[3]), []).append(i)

    for bare, idxs in by_bare.items():
        if bare not in PCG_BARE or len(idxs) != 1:
            continue  # only single-exon PCGs
        if bare in skip:
            continue  # owned by a reference-transfer pass (nad5-style join)
        i = idxs[0]
        r = rows[i]
        try:
            raw_s0, raw_e1 = int(r[1]), int(r[2])
        except (ValueError, IndexError):
            continue
        strand = r[5] if len(r) > 5 else "+"
        s0, e1 = _linear_span(raw_s0, raw_e1, n)

        nt = cds_from_row(gseq, raw_s0, raw_e1, strand)
        info = classify_cds(nt, code)
        if info["start_ok"] and info["stop_ok"] and not info["internal_stops"]:
            continue  # already a clean ORF

        def broken(reason):
            actions.append(f"PCG-skip {bare} ({reason})")
            unrepaired.append(bare)

        # Neighbour bounds in doubled coordinates, relaxed by the canonical
        # overlap this gene is allowed to share with its neighbour.
        slack = CANONICAL_OVERLAP.get(bare, 0)
        left_bound = max([sp[1] - slack for sp in spans if sp[1] <= s0] + [0])
        right_bound = min([sp[0] + slack for sp in spans if sp[0] >= e1] + [len(gseq2)])
        win_lo = max(left_bound, s0 - window)
        win_hi = min(right_bound, e1 + window)
        if win_lo >= win_hi:
            broken("no window")
            continue

        win = gseq2[win_lo:win_hi]
        if strand == "-":
            win = revcomp(win)
            s_from = (win_hi - e1) + 1
        else:
            s_from = (s0 - win_lo) + 1

        orf = refine_orf(win, s_from, len(win), code)
        if orf is None:
            broken("no clean ORF in window")
            continue
        cds_lo, cds_hi, aa, start_codon, _poly_a = orf
        old_len = e1 - s0
        new_len = cds_hi - cds_lo + 1
        if start_codon == "partial":
            broken(f"ORF not clean: start={start_codon}")
            continue
        floor = LEN_TOL_FLOOR.get(bare, DEFAULT_LEN_TOL_FLOOR)
        if abs(new_len - old_len) > max(floor, int(len_tol * old_len)):
            broken(f"len {old_len}->{new_len} outside tol")
            continue

        # Map window coords back to doubled-genome 0-based half-open, then onto
        # the real genome, restoring the start > end form for a wrapped CDS.
        if strand == "-":
            new_s0 = win_hi - cds_hi
            new_e1 = win_hi - cds_lo + 1
        else:
            new_s0 = win_lo + cds_lo - 1
            new_e1 = win_lo + cds_hi

        if new_s0 >= n:
            new_s0 -= n
            new_e1 -= n
        if new_e1 > n:
            new_e1 -= n          # wraps the origin: BED start > end, as MITOS writes it

        if (new_s0, new_e1) == (raw_s0, raw_e1):
            # The boundaries are already where refine_orf would put them, yet the
            # CDS did not classify clean -- an internal defect the snap can't fix.
            unrepaired.append(bare)
            actions.append(f"PCG-skip {bare} (no better boundary)")
            continue
        fixes[i] = (new_s0, new_e1)
        actions.append(
            f"PCG-fix {bare} {raw_s0}-{raw_e1}->{new_s0}-{new_e1} "
            f"({new_len}nt {len(aa)}aa)")

    return fixes, actions, unrepaired


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
    ap.add_argument("--pcg-window", type=int, default=36,
                    help="Half-width (nt) of the ORF-snap search window around a "
                         "mis-called PCG boundary, clipped at the neighbouring rows.")
    ap.add_argument("--pcg-len-tol", type=float, default=0.15,
                    help="A snapped PCG ORF whose length differs from the original "
                         "by more than max(60 nt, this fraction) is not applied.")
    ap.add_argument("--cox1-min-gain", type=int, default=30,
                    help="Replace MITOS's cox1 with the rebuilt intron-split join "
                         "when the clean ORF is at least this many nt longer. A "
                         "mis-called cox1 (the intron's homing-endonuclease ORF) "
                         "is typically barely half the real gene.")
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
        rrnl_seq, nad5_exons, cox1_seq = ref_features(args.ref_gb)
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
    # Repairs that were needed but could not be made. A sample with any of these
    # is PARTIAL, never FIXED: `add_rrnl is None` on its own cannot tell "16S was
    # already there" from "the transfer was rejected below thresholds", and
    # reporting the latter as FIXED is how a still-broken annotation gets waved
    # through downstream.
    rejected = []

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
            rejected.append("16S")
            actions.append(f"16S-skip cov={cov:.2f} pid={pid:.1f}")
    elif not have_rrnl:
        # 16S is missing and the reference probe did not even land.
        rejected.append("16S")
        actions.append("16S-skip no reference hit")

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
                    rejected.append("nad5")
                    actions.append(f"nad5-skip internal_stops={internal} (kept MITOS nad5)")
                else:
                    actions.append(f"nad5-ok (MITOS {cur_nt}nt; "
                                   f"rebuilt {ln if ok else 'NA'}nt clean={ok})")
            else:
                rejected.append("nad5")
                actions.append("nad5-skip below thresholds")

    # ---- cox1 group-I intron ----
    # Only worth attempting when MITOS's cox1 looks too short to be cox1: a
    # correctly called single-exon cox1 must not be second-guessed by a probe
    # transfer, and the reference's own cox1 is the yardstick for "too short".
    # An intron-split cox1 comes back at roughly half length, so this is a wide
    # margin, not a hair trigger.
    add_cox1 = None
    cox1_rows = [r for r in rows if bare_name(r[3]) == "cox1"]
    if cox1_seq and len(cox1_rows) == 1:
        r = cox1_rows[0]
        cur_nt = len(cds_from_row(gseq, int(r[1]), int(r[2]), r[5]))
        if cur_nt < len(cox1_seq) - args.cox1_min_gain:
            exons, nt, note = cox1_intron_join(
                cox1_seq, args.genome, (gseq + gseq), len(gseq), args.code,
                min_cov=args.min_cov, min_pid=args.min_pid)
            if exons is None:
                rejected.append("cox1")
                actions.append(f"cox1-skip {note} (kept MITOS cox1)")
            elif len(nt) - cur_nt < args.cox1_min_gain:
                rejected.append("cox1")
                actions.append(f"cox1-skip rebuilt {len(nt)}nt no better than "
                               f"MITOS {cur_nt}nt")
            else:
                add_cox1 = []
                spans = []
                for i, (s1, e1) in enumerate(exons):
                    # Map back off the doubled genome, restoring the start > end
                    # form MITOS uses for an origin-spanning feature.
                    b0 = (s1 - 1) % len(gseq)
                    b1 = e1 if e1 <= len(gseq) else e1 - len(gseq)
                    add_cox1.append([chrom, b0, b1, f"cox1_{i}", "0.0", "+"])
                    spans.append(f"{b0 + 1}-{b1}")
                actions.append(f"cox1 join {'+'.join(spans)} {len(nt) // 3 - 1}aa "
                               f"stops=0 (intron-split; MITOS had {cur_nt}nt, "
                               f"one exon) {note}")

    # ---- per-PCG ORF-boundary snap (assembly-only, no reference) ----
    # cox1 is excluded once the join owns it, exactly as nad5 always is.
    pcg_fixes, pcg_actions, pcg_unrepaired = snap_pcg_rows(
        rows, gseq, args.code, window=args.pcg_window, len_tol=args.pcg_len_tol,
        skip={"cox1"} if add_cox1 is not None else ())
    actions.extend(pcg_actions)
    rejected.extend(pcg_unrepaired)

    if add_rrnl is None and add_nad5 is None and add_cox1 is None and not pcg_fixes:
        # Nothing was applied. Still record what was found broken and left that
        # way, so SKIPPED is distinguishable from "there was nothing to do".
        if rejected:
            actions.append("unrepaired=" + ",".join(sorted(set(rejected))))
        finish("SKIPPED", ("; ".join(actions)) or "nothing to add", rows)

    # ---- build patched bed ----
    out_rows = []
    for i, r in enumerate(rows):
        bn = bare_name(r[3])
        if add_nad5 is not None and bn == "nad5":
            continue  # replaced below
        if add_cox1 is not None and bn == "cox1":
            continue  # MITOS's single-exon cox1; replaced by the join below
        if add_rrnl is not None and new_rrnl_span and TRNA_RE.match(r[3]):
            s0, e1 = int(r[1]) + 1, int(r[2])
            if s0 >= new_rrnl_span[0] and e1 <= new_rrnl_span[1]:
                actions.append(f"-{r[3]} (spurious tRNA inside 16S)")
                continue
        if i in pcg_fixes:
            new_s0, new_e1 = pcg_fixes[i]
            r = list(r)
            r[1], r[2] = str(new_s0), str(new_e1)
        out_rows.append(r)
    if add_rrnl is not None:
        out_rows.append(add_rrnl)
    if add_nad5 is not None:
        out_rows.extend(add_nad5)
    if add_cox1 is not None:
        out_rows.extend(add_cox1)

    # FIXED means "nothing this fixer was asked to repair is still broken".
    # Anything in `rejected` -- a 16S/nad5 transfer below thresholds, or a PCG
    # whose ORF the snap could not clean (OG2361's 873 nt cox1, OG2377's atp8
    # with no stop) -- keeps the sample PARTIAL so downstream still treats it as
    # deficient.
    if rejected:
        actions.append("unrepaired=" + ",".join(sorted(set(rejected))))
    state = "PARTIAL" if rejected else "FIXED"
    finish(state, "; ".join(actions), out_rows)


if __name__ == "__main__":
    main()
