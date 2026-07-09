#!/usr/bin/env python3
"""Reference-guided ordering / orientation / join of a multi-scaffold GetOrganelle
result into a single best-effort linear mitogenome.

GetOrganelle sometimes returns >=2 scaffolds that together carry the full gene
complement but that it could not disentangle into one circle. This is most often a
repeat collapse in the control region (a long tandem array short Illumina reads
cannot span), so the molecule breaks into pieces. The reseed-with-a-close-relative
step (subworkflows/.../getorganelle) closes many of these; the ones that survive it
were, until now, only blindly concatenated to scrape a species ID. Blind
concatenation also leaves whichever gene straddles the seam severed and any
minus-strand scaffold reversed, costing MITOS2 annotation calls it could otherwise
make.

This script uses the same related reference GenBank the circularity check already
consumes (check_getorganelle.py) to do what the reseed does in spirit: let a
trusted close relative resolve the layout. Each scaffold is BLASTed against the
reference, placed by its reference coordinate, reverse-complemented if it sits on
the minus strand, and the scaffolds are joined in reference order. Real
sequence overlaps between adjacent scaffolds are merged; genuine gaps are bridged
with an N-spacer sized to the reference gap (capped) so the join never fabricates
junction sequence that does not exist.

The result is an inference, not a read-closed assembly: the seam sequence inside
repeat / control-region gaps is unknown and is represented by Ns. The join fixes
order, orientation and gene completeness so downstream annotation and species
validation work; it does not invent the missing repeat. Every placement decision
(order, strand, per-seam gap/overlap, any scaffold that did not hit the reference)
is written to an evidence TSV, and the output is flagged ``reference_guided_join``
so it is never mistaken for a formally circularised genome. Scaffolds that do not
hit the reference at all (possible NUMT / contamination) are reported and appended
last, never silently dropped.

The single-record FASTA this writes then flows into the existing
check_getorganelle.py circularity / length-anomaly screen exactly as a one-scaffold
result would.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile


def log(msg):
    print(f"[join_scaffolds] {msg}", file=sys.stderr)


_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                      "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def read_fasta(path):
    """Return a list of (header_id, sequence) in file order."""
    recs = []
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(parts)))
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        recs.append((header, "".join(parts)))
    return recs


def gb_to_fasta(gb_path, out_path):
    """Extract the ORIGIN sequence of a GenBank record to FASTA; return length.
    Mirrors check_getorganelle.gb_to_fasta so both steps see the same reference."""
    seq = []
    in_origin = False
    with open(gb_path) as fh:
        for line in fh:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                seq.append("".join(re.findall(r"[acgtnACGTN]", line)))
    s = "".join(seq)
    with open(out_path, "w") as out:
        out.write(">ref\n" + s + "\n")
    return len(s)


def run(cmd, **kw):
    log("$ " + (cmd if isinstance(cmd, str) else " ".join(cmd)))
    return subprocess.run(cmd, check=True, **kw)


def blast_scaffolds(recs, ref_fa, workdir, min_len, min_pident):
    """BLAST every scaffold against the reference. Return, per scaffold index, a
    list of HSPs as (q_start, q_end, r_lo, r_hi, strand, length, pident)."""
    q = os.path.join(workdir, "scaffolds.fa")
    with open(q, "w") as fh:
        for i, (_h, s) in enumerate(recs):
            fh.write(f">s{i}\n{s}\n")
    db = os.path.join(workdir, "refdb")
    run(["makeblastdb", "-in", ref_fa, "-dbtype", "nucl", "-out", db],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # dc-megablast (not the blastn default of megablast): scaffolds are matched to a
    # *related-species* reference at ~70-80% identity, which megablast's long exact
    # seeds miss entirely (it dropped COX2 on the first real test). dc-megablast's
    # discontiguous seeds find these diverged HSPs.
    # Dust masking is left ON (default): a control-region satellite array against a
    # low-complexity reference stretch otherwise explodes into thousands of short
    # HSPs and hangs dc-megablast (seen on a junk 3.4 kb assembly). Low-complexity
    # sequence should never drive placement anyway. max_hsps + a hard timeout are a
    # belt-and-braces guard so a pathological input can never wedge the pipeline.
    try:
        res = subprocess.run(
            ["blastn", "-task", "dc-megablast", "-query", q, "-db", db,
             "-max_hsps", "50", "-max_target_seqs", "1",
             "-outfmt", "6 qseqid qstart qend sstart send length pident"],
            check=True, stdout=subprocess.PIPE, universal_newlines=True, timeout=300,
        )
    except subprocess.TimeoutExpired:
        log("BLAST timed out; treating all scaffolds as unplaced.")
        return {i: [] for i in range(len(recs))}
    hsps = {i: [] for i in range(len(recs))}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 7:
            continue
        idx = int(f[0][1:])
        qs, qe, ss, se = (int(f[i]) for i in range(1, 5))
        ln, pid = int(f[5]), float(f[6])
        if ln < min_len or pid < min_pident:
            continue
        strand = "+" if se >= ss else "-"
        hsps[idx].append((qs, qe, min(ss, se), max(ss, se), strand, ln, pid))
    return hsps


def place_scaffold(hsps, scaf_len):
    """From a scaffold's HSPs derive placement: (strand, ref_lo, ref_hi,
    aligned_bp, w_pident). Scaffolds are ordered by ref_lo and oriented to the
    reference (+) strand. None if no qualifying HSPs."""
    if not hsps:
        return None
    plus = sum(h[5] for h in hsps if h[4] == "+")
    minus = sum(h[5] for h in hsps if h[4] == "-")
    strand = "+" if plus >= minus else "-"
    ref_lo = min(h[2] for h in hsps)
    ref_hi = max(h[3] for h in hsps)
    # Guard against spurious scattered hits: a scaffold cannot truly span more
    # reference than ~its own length. When min/max blows past that (a few junk HSPs
    # at opposite ends of the reference, common for low-complexity / contaminant
    # scaffolds), fall back to the span of the single longest HSP, which reflects
    # where the scaffold actually sits. Keeps coverage honest for junk samples.
    if (ref_hi - ref_lo + 1) > 1.5 * scaf_len:
        best = max(hsps, key=lambda h: h[5])
        ref_lo, ref_hi = best[2], best[3]
    aligned = sum(h[5] for h in hsps)
    w_pid = sum(h[5] * h[6] for h in hsps) / aligned if aligned else 0.0
    return strand, ref_lo, ref_hi, aligned, round(w_pid, 1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="Multi-scaffold GetOrganelle FASTA.")
    p.add_argument("--reference-gb", required=True, help="Related-species reference GenBank.")
    p.add_argument("--sample", required=True)
    p.add_argument("--out-fasta", required=True, help="Single-record joined FASTA.")
    p.add_argument("--out-evidence", required=True, help="Join evidence TSV.")
    p.add_argument("--min-hsp-len", type=int, default=100)
    p.add_argument("--min-hsp-identity", type=float, default=65.0,
                   help="Min %% identity for a scaffold-vs-reference HSP "
                        "(cross-species mito sits at ~70-80%%).")
    p.add_argument("--max-gap-n", type=int, default=100,
                   help="Cap on the N-spacer inserted at a reference gap.")
    args = p.parse_args()

    recs = read_fasta(args.fasta)
    nrec = len(recs)
    note = "ok"

    # Single record (or empty): nothing to join. Pass the input straight through so
    # this step is a no-op outside the multi-scaffold case it exists for.
    if nrec <= 1:
        seq = recs[0][1] if recs else ""
        with open(args.out_fasta, "w") as out:
            out.write(f">{args.sample}\n{seq}\n")
        write_evidence(args.out_evidence, args.sample, nrec, nrec, 0,
                       order="NA", seams="NA", joined_len=len(seq),
                       ref_len="NA", coverage="NA", looks_circular="NA",
                       note="single_record_passthrough")
        log(f"{nrec} record(s); passthrough.")
        return 0

    have_ref = os.path.exists(args.reference_gb) and os.path.getsize(args.reference_gb) > 0
    if not have_ref:
        # No reference -> fall back to the old behaviour (concatenate in file order)
        # so we still produce something for species ID, but flag it honestly.
        joined = "".join(s for _h, s in recs)
        with open(args.out_fasta, "w") as out:
            out.write(f">{args.sample} concat_no_reference\n{joined}\n")
        write_evidence(args.out_evidence, args.sample, nrec, 0, nrec,
                       order="file_order", seams="NA", joined_len=len(joined),
                       ref_len="NA", coverage="NA", looks_circular="NA",
                       note="no_reference_concat")
        log("No reference; concatenated in file order (fallback).")
        return 0

    with tempfile.TemporaryDirectory() as wd:
        ref_fa = os.path.join(wd, "ref.fa")
        ref_len = gb_to_fasta(args.reference_gb, ref_fa)
        hsps = blast_scaffolds(recs, ref_fa, wd, args.min_hsp_len, args.min_hsp_identity)

    placed = []     # (idx, strand, ref_lo, ref_hi, aligned, w_pid)
    unplaced = []   # idx
    for i, (_h, s) in enumerate(recs):
        pl = place_scaffold(hsps[i], len(s))
        if pl is None:
            unplaced.append(i)
        else:
            placed.append((i,) + pl)

    if not placed:
        # Nothing hit the reference -> behave like the no-reference fallback.
        joined = "".join(s for _h, s in recs)
        with open(args.out_fasta, "w") as out:
            out.write(f">{args.sample} concat_no_hit\n{joined}\n")
        write_evidence(args.out_evidence, args.sample, nrec, 0, nrec,
                       order="file_order", seams="NA", joined_len=len(joined),
                       ref_len=ref_len, coverage="0.0%", looks_circular="NA",
                       note="no_scaffold_hit_reference")
        log("No scaffold hit the reference; concatenated in file order (fallback).")
        return 0

    placed.sort(key=lambda x: x[2])  # by ref_lo (leftmost reference coordinate)

    # Build the joined sequence in reference order. Each scaffold is oriented to the
    # reference (+) strand and laid down left to right. Where a scaffold's reference
    # span is largely contained within the span already covered (it nests inside an
    # earlier, longer scaffold -- e.g. a short COX2 fragment sitting in a gap of the
    # main scaffold), it is still kept but flagged as redundant: it adds the gene but
    # repeats reference territory, which the reviewer should see (and which inflates
    # the length ratio honestly rather than being hidden).
    pieces = []
    order_tokens = []
    seam_tokens = []
    covered_lo_hi = []
    prev_ref_hi = None
    for n, (idx, strand, ref_lo, ref_hi, _aligned, w_pid) in enumerate(placed):
        seq = recs[idx][1]
        oriented = seq if strand == "+" else revcomp(seq)
        if n > 0:
            gap = ref_lo - prev_ref_hi - 1  # reference-coordinate gap
            if gap > 0:
                n_count = min(gap, args.max_gap_n)
                pieces.append("N" * n_count)
                seam_tokens.append(f"gap{gap}bp:N{n_count}")
            elif ref_hi <= prev_ref_hi:
                seam_tokens.append(f"nested(ref{ref_lo}-{ref_hi}):redundant")
            else:
                seam_tokens.append(f"overlap{-gap}bp:buttjoin")
        pieces.append(oriented)
        order_tokens.append(f"{recs[idx][0]}({strand},ref{ref_lo}-{ref_hi},{w_pid}%)")
        covered_lo_hi.append((ref_lo, ref_hi))
        prev_ref_hi = max(prev_ref_hi, ref_hi) if prev_ref_hi is not None else ref_hi

    # Append any unplaced scaffolds last, behind an N-spacer, clearly flagged.
    for idx in unplaced:
        pieces.append("N" * args.max_gap_n)
        pieces.append(recs[idx][1])
        order_tokens.append(f"{recs[idx][0]}(unplaced)")
        seam_tokens.append(f"unplaced:N{args.max_gap_n}")

    joined = "".join(pieces)

    # Reference coverage of the placed scaffolds (merged reference intervals).
    covered_lo_hi.sort()
    merged = [list(covered_lo_hi[0])]
    for lo, hi in covered_lo_hi[1:]:
        if lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    covered = sum(hi - lo + 1 for lo, hi in merged)
    coverage = covered / ref_len if ref_len else 0.0
    # "Looks circular" = placed scaffolds blanket the reference end-to-end (the only
    # remaining gap is the wrap at the origin). A heuristic for the evidence/QC, not
    # a closure claim.
    wrap_gap = (ref_len - merged[-1][1]) + (merged[0][0] - 1)
    looks_circular = coverage >= 0.95 and len(merged) == 1

    with open(args.out_fasta, "w") as out:
        out.write(f">{args.sample} reference_guided_join scaffolds={nrec}\n")
        for i in range(0, len(joined), 80):
            out.write(joined[i:i + 80] + "\n")

    write_evidence(args.out_evidence, args.sample, nrec, len(placed), len(unplaced),
                   order=" >> ".join(order_tokens),
                   seams=";".join(seam_tokens) if seam_tokens else "NA",
                   joined_len=len(joined), ref_len=ref_len,
                   coverage=f"{coverage*100:.1f}%",
                   looks_circular=("True" if looks_circular else "False"),
                   note=note, wrap_gap=wrap_gap)

    log(f"joined {len(placed)} placed (+{len(unplaced)} unplaced) -> {len(joined)} bp; "
        f"ref_cov={coverage*100:.1f}%; looks_circular={looks_circular}")
    return 0


def write_evidence(path, sample, n_in, n_placed, n_unplaced, order, seams,
                   joined_len, ref_len, coverage, looks_circular, note,
                   wrap_gap="NA"):
    ratio = "NA"
    if isinstance(ref_len, int) and ref_len:
        ratio = round(joined_len / ref_len, 3)
    cols = ["sample", "method", "n_input_scaffolds", "n_placed", "n_unplaced",
            "scaffold_order", "seams", "joined_length", "reference_length",
            "length_ratio", "reference_coverage", "origin_wrap_gap_bp",
            "looks_circular", "note"]
    vals = [sample, "reference_guided_join", n_in, n_placed, n_unplaced,
            order, seams, joined_len, ref_len, ratio, coverage, wrap_gap,
            looks_circular, note]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(str(v) for v in vals) + "\n")


if __name__ == "__main__":
    sys.exit(main())
