#!/usr/bin/env python3
"""Collapse a clean tandem-genome concatemer to a single monomer.

MitoHiFi / GetOrganelle occasionally emit a head-to-tail dimer (or multimer):
a circular mitogenome assembled at ~2x (or 3x) the true length. The circularity
check (bin/check_circularity.py / bin/check_getorganelle.py) already detects this
and records anomaly_type=concatemer plus a suggested_trim_region; this script
*applies* that curation so the assembly finishes automatically instead of going
to manual review.

Only a genuine, near-integer multimer is collapsed. The suggested monomer is
re-verified by self-alignment (the monomer must cover the trimmed tail at high
identity) before anything is rewritten; anything that fails the check is passed
through unchanged with an explanatory report, so the step can never silently
corrupt an assembly it does not understand.
"""
import argparse
import os
import subprocess
import sys
import tempfile


def log(msg):
    print(f"[collapse_concatemer] {msg}", file=sys.stderr)


def read_fasta(path):
    """Return (header_line_without_gt, sequence) for the first record."""
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                else:
                    break
            elif header is not None:
                seq.append(line.strip())
    return header, "".join(seq)


def read_evidence(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if len(lines) < 2:
        return {}
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    return {h.strip(): (values[i].strip() if i < len(values) else "") for i, h in enumerate(header)}


def write_corrected_evidence(source_path, out_path, action, collapsed_length, period, coverage):
    """Preserve the original check evidence and append the curation outcome.

    Downstream QC consumes this corrected sidecar.  The original check remains
    published unchanged for provenance.
    """
    with open(source_path) as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    header = lines[0].split("\t") if lines else []
    values = lines[1].split("\t") if len(lines) > 1 else []
    while len(values) < len(header):
        values.append("")
    row = dict(zip(header, values))
    original_anomaly = row.get("anomaly_type", "")
    original_breakpoint = row.get("suggested_trim_region", "")
    original_repeat_identity = row.get("repeat_identity", "")
    if action == "collapsed":
        row["assembly_length"] = str(collapsed_length)
        row["length_ratio"] = (
            str(round(collapsed_length / int(row["reference_length"]), 3))
            if (row.get("reference_length") or "").isdigit() else "NA"
        )
        row["excess_bp"] = "0"
        row["length_anomaly"] = "no"
        row["anomaly_type"] = "none"
        row["suggested_trim_region"] = "NA"
        row["curation_suggestion"] = "resolved_by_concatemer_collapse"
    extra = {
        "curation_action": action,
        "pre_curation_anomaly_type": original_anomaly or "none",
        "pre_curation_suggested_trim_region": original_breakpoint or "NA",
        "pre_curation_repeat_identity": original_repeat_identity or "NA",
        "inferred_monomer_period": str(period) if period else "NA",
        "curation_verification_coverage": (
            str(round(coverage, 3)) if coverage is not None else "NA"
        ),
    }
    for key, value in extra.items():
        if key not in header:
            header.append(key)
        row[key] = value
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        out.write("\t".join(str(row.get(col, "")) for col in header) + "\n")


def parse_trim(region):
    """'start-end' (1-based, inclusive) -> (start, end) or None."""
    try:
        start, end = region.split("-")
        return int(start), int(end)
    except (ValueError, AttributeError):
        return None


def self_identity(monomer, tail, workdir, expected_length=None):
    """Fraction of `tail` covered by `monomer` under blastn. A true concatemer
    tail is ~identical to the monomer; a chimeric over-length assembly is not."""
    if not monomer or not tail:
        return 0.0
    qf = os.path.join(workdir, "tail.fa")
    df = os.path.join(workdir, "mono.fa")
    with open(qf, "w") as fh:
        fh.write(">tail\n" + tail + "\n")
    with open(df, "w") as fh:
        fh.write(">mono\n" + monomer + "\n")
    db = os.path.join(workdir, "monodb")
    subprocess.run(["makeblastdb", "-in", df, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res = subprocess.run(
        ["blastn", "-task", "megablast", "-query", qf, "-db", db,
         "-outfmt", "6 qstart qend pident"],
        check=True, stdout=subprocess.PIPE, universal_newlines=True,
    )
    covered = [False] * len(tail)
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and float(f[2]) >= 90.0:
            a, b = int(f[0]), int(f[1])
            for i in range(min(a, b) - 1, max(a, b)):
                if 0 <= i < len(tail):
                    covered[i] = True
    denominator = min(len(tail), expected_length or len(tail))
    return min(1.0, sum(covered) / max(denominator, 1))


def infer_monomer_period(seq, ref_len, workdir):
    """Infer a tandem-genome period from long off-diagonal self alignments.

    The reference length is a safety window, not the breakpoint.  This matters
    for assemblies such as OG750 where the two assembled copies are ~16.3 kb but
    the closest reference is 16.7 kb.
    """
    if not seq or not ref_len:
        return [], "missing_reference_length"
    seq_file = os.path.join(workdir, "assembly.fa")
    with open(seq_file, "w") as fh:
        fh.write(">assembly\n" + seq + "\n")
    db = os.path.join(workdir, "assemblydb")
    subprocess.run(
        ["makeblastdb", "-in", seq_file, "-dbtype", "nucl", "-out", db],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    result = subprocess.run(
        ["blastn", "-task", "megablast", "-dust", "no", "-query", seq_file,
         "-db", db, "-outfmt", "6 qstart qend sstart send length pident"],
        check=True, stdout=subprocess.PIPE, universal_newlines=True,
    )
    lower, upper = 0.9 * ref_len, 1.15 * ref_len
    scores = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        qs, qe, ss, se, aln_len = map(int, fields[:5])
        pident = float(fields[5])
        if pident < 90.0 or (qe - qs) * (se - ss) <= 0:
            continue
        period = abs(ss - qs)
        if period == 0 or not (lower <= period <= upper) or aln_len < 0.25 * ref_len:
            continue
        # Long genome-scale HSPs dominate short control-region repeats.
        scores[period] = scores.get(period, 0.0) + aln_len * (pident / 100.0)
    if not scores:
        return [], "no_reference_scale_off_diagonal_alignment"
    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    return [period for period, _score in ranked], "ok"


def write_fasta(path, header, seq, wrap=60):
    with open(path, "w") as out:
        out.write(">" + header + "\n")
        for i in range(0, len(seq), wrap):
            out.write(seq[i:i + wrap] + "\n")


def write_report(path, fields):
    cols = ["sample", "action", "original_length", "collapsed_length",
            "reference_length", "inferred_monomer_period", "tail_identity", "reason"]
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(str(fields.get(c, "NA")) for c in cols) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--evidence", required=True, help="circularity_check.tsv from the check step")
    p.add_argument("--sample", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-report", required=True)
    p.add_argument("--out-evidence", required=True,
                   help="Corrected post-curation evidence consumed by downstream QC.")
    p.add_argument("--min-tail-identity", type=float, default=0.9,
                   help="Fraction of the trimmed tail the monomer must cover to confirm a concatemer.")
    args = p.parse_args()

    header, seq = read_fasta(args.fasta)
    ev = read_evidence(args.evidence)
    anomaly = (ev.get("anomaly_type") or "").strip().lower()
    ref_len = ev.get("reference_length")
    ref_len = int(ref_len) if (ref_len or "").isdigit() else None
    report = {"sample": args.sample, "action": "passthrough",
              "original_length": len(seq), "collapsed_length": len(seq),
              "reference_length": ref_len if ref_len else "NA",
              "inferred_monomer_period": "NA", "tail_identity": "NA", "reason": ""}

    if header is None or not seq:
        report["reason"] = "no_sequence"
        write_fasta(args.out_fasta, header or args.sample, seq)
        write_report(args.out_report, report)
        write_corrected_evidence(args.evidence, args.out_evidence, "passthrough", len(seq), None, None)
        return 0

    if anomaly != "concatemer":
        report["reason"] = f"not_a_concatemer(anomaly={anomaly or 'none'})"
        write_fasta(args.out_fasta, header, seq)
        write_report(args.out_report, report)
        write_corrected_evidence(args.evidence, args.out_evidence, "passthrough", len(seq), None, None)
        log(report["reason"])
        return 0

    try:
        with tempfile.TemporaryDirectory() as wd:
            periods, period_reason = infer_monomer_period(seq, ref_len, wd)
            verified = []
            for candidate in periods[:10]:
                coverage = self_identity(seq[:candidate], seq[candidate:], wd, ref_len)
                if coverage >= args.min_tail_identity:
                    verified.append((candidate, coverage))
    except subprocess.CalledProcessError as e:
        verified, period_reason = [], f"period_blast_failed({e})"
    if not verified:
        report["reason"] = period_reason
        write_fasta(args.out_fasta, header, seq)
        write_report(args.out_report, report)
        write_corrected_evidence(args.evidence, args.out_evidence, "passthrough", len(seq), None, None)
        log(report["reason"])
        return 0

    verified.sort(key=lambda item: (-item[1], item[0]))
    period, ident = verified[0]
    if any(abs(other_period - period) > 50 and other_coverage >= ident - 0.02
           for other_period, other_coverage in verified[1:]):
        report["reason"] = "ambiguous_monomer_period"
        write_fasta(args.out_fasta, header, seq)
        write_report(args.out_report, report)
        write_corrected_evidence(args.evidence, args.out_evidence, "passthrough", len(seq), None, None)
        log(report["reason"])
        return 0

    report["inferred_monomer_period"] = period
    monomer = seq[:period]
    tail = seq[period:]

    # Re-verify: monomer length near the reference, and the tail is (near-)identical
    # to the monomer. Fail closed to passthrough if either check does not hold.
    if ref_len and not (0.9 * ref_len <= len(monomer) <= 1.15 * ref_len):
        report["reason"] = f"monomer_length_off_reference({len(monomer)}vs{ref_len})"
        write_fasta(args.out_fasta, header, seq)
        write_report(args.out_report, report)
        write_corrected_evidence(args.evidence, args.out_evidence, "passthrough", len(seq), period, None)
        log(report["reason"])
        return 0

    report["tail_identity"] = round(ident, 3)
    if ident < args.min_tail_identity:
        report["reason"] = f"tail_identity_below_threshold({round(ident,3)})"
        write_fasta(args.out_fasta, header, seq)
        write_report(args.out_report, report)
        write_corrected_evidence(args.evidence, args.out_evidence, "passthrough", len(seq), period, ident)
        log(report["reason"])
        return 0

    report["action"] = "collapsed"
    report["collapsed_length"] = len(monomer)
    report["reason"] = f"collapsed_{round(len(seq)/max(len(monomer),1),2)}x_concatemer"
    write_fasta(args.out_fasta, header, monomer)
    write_report(args.out_report, report)
    write_corrected_evidence(args.evidence, args.out_evidence, "collapsed", len(monomer), period, ident)
    log(f"collapsed {len(seq)} -> {len(monomer)} bp (tail identity {round(ident,3)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
