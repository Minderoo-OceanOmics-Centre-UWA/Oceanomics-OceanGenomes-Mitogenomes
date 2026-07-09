#!/usr/bin/env python3
"""Re-test circularity and screen length/tandem-repeat anomalies for a
GetOrganelle assembly, mirroring the MitoHiFi check (bin/check_circularity.py).

GetOrganelle reports a result status of "circular genome" or "N scaffold(s)".
A single non-circularised scaffold is frequently still a *complete* circle that
GetOrganelle just could not formally close (its docs note this). We re-test such
scaffolds against the related complete reference: a complete circle linearised at
a different origin covers ~all of the reference, whereas a genuinely incomplete
scaffold leaves a contiguous chunk of the reference uncovered. When the reference
test confirms a full circle, the circular verdict is corrected to True; the
original GetOrganelle status and the supporting evidence are kept in the
evidence TSV.

The same run also screens length / tandem-repeat anomalies (concatemer,
control-region repeat, unresolved over-length) exactly as the HiFi check does, so
GetOrganelle assemblies feed the same QC-gate anomaly block and assembly-summary
manual-review reason. CR location needs the assembly's own annotation, which is
not available at assembly time, so repeat_in_control_region is left NA here; the
concatemer / length-outlier classification does not depend on it.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile


def log(msg):
    print(f"[check_getorganelle] {msg}", file=sys.stderr)


def read_records(path):
    """Return (num_records, first_record_sequence, total_length)."""
    nrec = 0
    first = []
    total = 0
    cur_is_first = False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                nrec += 1
                cur_is_first = nrec == 1
            else:
                s = line.strip()
                total += len(s)
                if cur_is_first:
                    first.append(s)
    return nrec, "".join(first), total


def gb_to_fasta(gb_path, out_path):
    """Extract the ORIGIN sequence of a GenBank record to FASTA; return length."""
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


def blast_reference_coverage(seq, ref_fa, ref_len, workdir):
    """BLAST the assembly against the reference; return (coverage_fraction,
    biggest_uncovered_gap_bp)."""
    q = os.path.join(workdir, "asm.fa")
    with open(q, "w") as fh:
        fh.write(">asm\n" + seq + "\n")
    db = os.path.join(workdir, "refdb")
    run(["makeblastdb", "-in", ref_fa, "-dbtype", "nucl", "-out", db],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # dc-megablast, not the blastn default of megablast: the reference is a *related
    # species* at ~70-80% identity, and megablast's long exact seeds miss most of
    # those diverged HSPs -- it reported 1.3% reference coverage on a scaffold that
    # dc-megablast covers at 93%, which would wrongly leave a complete circle uncalled.
    res = subprocess.run(
        ["blastn", "-task", "dc-megablast", "-query", q, "-db", db,
         "-outfmt", "6 sstart send length pident"],
        check=True, stdout=subprocess.PIPE, universal_newlines=True,
    )
    hsps = []
    for line in res.stdout.splitlines():
        f = line.split("\t")
        # 70% identity floor: a related-species reference aligns at ~70-80%, so an
        # 80% floor (with dc-megablast) would still discard most genuine coverage.
        if len(f) >= 4 and int(f[2]) >= 100 and float(f[3]) >= 70:
            hsps.append((min(int(f[0]), int(f[1])), max(int(f[0]), int(f[1]))))
    if not hsps:
        return 0.0, ref_len
    hsps.sort()
    merged = [list(hsps[0])]
    for s, e in hsps[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    covered = sum(e - s + 1 for s, e in merged)
    gaps = [(b[0] - a[1] - 1) for a, b in zip(merged, merged[1:])]
    wrap = (ref_len - merged[-1][1]) + (merged[0][0] - 1)
    biggest = max(gaps + [wrap]) if (gaps or wrap) else 0
    return covered / ref_len, biggest


_REF_CONSUMING = set("MDN=X")


def find_tandem_repeat(seq, workdir, min_len, min_pident):
    """Self-align the assembly to locate a forward tandem-repeat array.
    Mirrors check_circularity.find_tandem_repeat."""
    fa = os.path.join(workdir, "self.fa")
    with open(fa, "w") as fh:
        fh.write(">s\n" + seq + "\n")
    db = os.path.join(workdir, "selfdb")
    run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res = subprocess.run(
        ["blastn", "-query", fa, "-db", db, "-dust", "no",
         "-outfmt", "6 qstart qend sstart send length pident"],
        check=True, stdout=subprocess.PIPE, universal_newlines=True,
    )
    hits = []
    for line in res.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 6:
            continue
        qs, qe, ss, se = (int(p[i]) for i in range(4))
        ln, pid = int(p[4]), float(p[5])
        if ln < min_len or pid < min_pident:
            continue
        if se <= ss or ss <= qs:
            continue
        hits.append((qs, qe, ss, se, ln, pid))
    if not hits:
        return None
    start = min(min(h[0], h[2]) for h in hits)
    end = max(max(h[1], h[3]) for h in hits)
    unit = min(h[2] - h[0] for h in hits)
    array_len = end - start + 1
    if unit <= 0 or array_len < max(2 * unit, min_len):
        return None
    tw = sum(h[4] for h in hits)
    ident = sum(h[4] * h[5] for h in hits) / tw if tw else 0.0
    return {"array_start": start, "array_end": end, "array_len": array_len,
            "unit_bp": unit, "copies": round(array_len / unit, 1), "identity": round(ident, 1)}


def coerce_bool(v):
    s = str(v).strip().lower()
    if s in {"true", "t", "1", "yes", "y", "circular"}:
        return True
    if s in {"false", "f", "0", "no", "n", "linear", "scaffold"}:
        return False
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--reference-gb", help="Related-species reference GenBank (for length + circularity test).")
    p.add_argument("--sample", required=True)
    p.add_argument("--getorg-circular", default="null", help="GetOrganelle log verdict: true / false / null.")
    p.add_argument("--out-evidence", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-ref-coverage", type=float, default=0.95, help="Reference coverage to call a scaffold a complete circle.")
    p.add_argument("--max-length-ratio", type=float, default=1.10)
    p.add_argument("--min-repeat-length", type=int, default=200)
    p.add_argument("--min-repeat-identity", type=float, default=70.0)
    args = p.parse_args()

    getorg_circ = coerce_bool(args.getorg_circular)
    nrec, seq, total_len = read_records(args.fasta)
    L = len(seq)
    note = "ok"

    ref_len = None
    ref_cov = None
    ref_gap = None
    circular_by_reference = None

    have_ref = args.reference_gb and os.path.exists(args.reference_gb) and os.path.getsize(args.reference_gb) > 0
    if have_ref and L:
        try:
            with tempfile.TemporaryDirectory() as wd:
                ref_fa = os.path.join(wd, "ref.fa")
                ref_len = gb_to_fasta(args.reference_gb, ref_fa)
                if ref_len > 0:
                    ref_cov, ref_gap = blast_reference_coverage(seq, ref_fa, ref_len, wd)
                    ratio = L / ref_len
                    circular_by_reference = (
                        nrec == 1 and ref_cov >= args.min_ref_coverage and 0.9 <= ratio <= 1.15
                    )
        except subprocess.CalledProcessError as e:
            note = "blast_failed"
            log(f"reference BLAST failed ({e})")
        except Exception as e:  # pragma: no cover
            note = "reference_error"
            log(f"reference assessment error ({e})")
    elif not have_ref:
        note = "no_reference"

    # Final circular verdict: keep GetOrganelle's True; otherwise adopt the
    # reference-confirmed circle. None stays None (unknown).
    if getorg_circ is True:
        final_circular = True
    elif circular_by_reference is True:
        final_circular = True
    elif getorg_circ is False or circular_by_reference is False:
        final_circular = False
    else:
        final_circular = None
    circular_corrected = bool(getorg_circ is not True and final_circular is True)

    # Length / tandem-repeat anomaly screen.
    length_ratio = excess = length_anomaly = None
    if ref_len:
        length_ratio = round(L / ref_len, 3)
        excess = L - ref_len
        length_anomaly = (L / ref_len) >= args.max_length_ratio

    repeat = None
    if seq:
        try:
            with tempfile.TemporaryDirectory() as wd:
                repeat = find_tandem_repeat(seq, wd, args.min_repeat_length, args.min_repeat_identity)
        except Exception as e:
            log(f"repeat scan error ({e})")

    frac = (repeat["array_len"] / L) if (repeat and L) else 0.0
    anomaly_type = "none"
    repeat_region = "NA"
    suggested_trim = "NA"
    curation = "none"
    if repeat:
        repeat_region = f"{repeat['array_start']}-{repeat['array_end']}"
    if length_anomaly and repeat and frac >= 0.8 and ref_len and L / ref_len >= 1.6:
        anomaly_type = "concatemer"
        suggested_trim = f"{ref_len + 1}-{L}"
        curation = (f"Assembly is {round(L / ref_len, 2)}x the reference (tandem genome duplication / "
                    f"concatemer). Collapse to a single monomer (~{ref_len} bp), removing {suggested_trim}.")
    elif length_anomaly and repeat:
        anomaly_type = "tandem_repeat"
        trim_bp = min(excess, repeat["array_len"] - repeat["unit_bp"])
        if trim_bp > 0:
            suggested_trim = f"{repeat['array_end'] - trim_bp + 1}-{repeat['array_end']}"
        curation = (f"Assembly {excess} bp longer than reference; tandem repeat (~{repeat['unit_bp']} bp unit "
                    f"x ~{repeat['copies']} copies) at {repeat_region}. Review for length heteroplasmy.")
    elif length_anomaly:
        anomaly_type = "unresolved"
        curation = f"Assembly {excess} bp longer than reference with no clear tandem repeat; review for duplication / NUMT."

    def fmt(v):
        return "NA" if v is None else ("True" if v is True else ("False" if v is False else str(v)))

    cols = ["sample", "getorg_circular", "num_records", "reference_coverage", "reference_max_gap",
            "circular_by_reference", "final_verdict_circular", "circular_corrected",
            "assembly_length", "reference_length", "length_ratio", "excess_bp", "length_anomaly",
            "tandem_repeat", "repeat_region", "anomaly_type", "suggested_trim_region",
            "curation_suggestion", "note"]
    vals = [args.sample, fmt(getorg_circ), nrec,
            (f"{ref_cov*100:.1f}%" if ref_cov is not None else "NA"),
            (ref_gap if ref_gap is not None else "NA"),
            fmt(circular_by_reference), fmt(final_circular), ("yes" if circular_corrected else "no"),
            L, (ref_len if ref_len else "NA"), fmt(length_ratio), fmt(excess),
            ("yes" if length_anomaly else ("no" if length_anomaly is False else "NA")),
            ("yes" if repeat else "no"), repeat_region, anomaly_type, suggested_trim, curation, note]

    os.makedirs(os.path.dirname(args.out_evidence) or ".", exist_ok=True)
    with open(args.out_evidence, "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(str(v) for v in vals) + "\n")

    log(f"circular: getorg={fmt(getorg_circ)} ref={fmt(circular_by_reference)} -> final={fmt(final_circular)} "
        f"(corrected={circular_corrected}); anomaly={anomaly_type}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
