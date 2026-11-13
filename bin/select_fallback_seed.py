#!/usr/bin/env python3
"""Choose a bounded, read-supported seed panel when GetOrganelle's first pass is empty.

The normal invertebrate selector ranks a curated group database against the first-pass
assembly. An empty assembly has no sequence to rank, so this selector first narrows the
manifest at the most specific available taxonomic level and then ranks that shortlist
with a fixed read subsample. It always emits an audit status and exits zero.
"""
import argparse
import gzip
import math
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# bin/ is on PATH at runtime, so a sibling import resolves via sys.path[0].
# Bio is imported lazily inside the two functions that need it (same pattern as
# refdb_record) so the taxonomy shortlist stays unit-testable without Biopython.
import refdb_record
from refdb_record import subset_label_db, union_coverage

UNRESOLVED = {"", "unknown", "na", "none", "dropped"}


def norm(value):
    return re.sub(r"[^a-z0-9]+", " ", str(value or "").casefold()).strip()


def record_terms(row):
    return {norm(row.get("organism")), norm(row.get("family")),
            *(norm(x) for x in row.get("lineage", []))}


def balanced_cap(accessions, rows, limit):
    """Deterministically retain broad family representation within a large tier."""
    if len(accessions) <= limit:
        return sorted(accessions)
    by_family = defaultdict(list)
    for acc in sorted(accessions):
        by_family[norm(rows[acc].get("family")) or "_unknown"].append(acc)
    result = []
    while len(result) < limit and by_family:
        for family in sorted(list(by_family)):
            result.append(by_family[family].pop(0))
            if not by_family[family]:
                del by_family[family]
            if len(result) == limit:
                break
    return result


def taxonomy_shortlist(rows, taxon="", family="", order="", tax_class="", max_candidates=50):
    queries = [("nominal", taxon), ("family", family), ("order", order), ("class", tax_class)]
    missed = []
    for level, value in queries:
        key = norm(value)
        if key in UNRESOLVED:
            continue
        hits = [acc for acc, row in rows.items() if key in record_terms(row)]
        if hits:
            return balanced_cap(hits, rows, max_candidates), level, missed
        missed.append(f"{level}:{value}")
    return balanced_cap(list(rows), rows, max_candidates), "group", missed


def open_reads(path):
    handle = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    suffixes = Path(str(path).removesuffix(".gz")).suffixes
    fmt = "fastq" if any(s in {".fq", ".fastq"} for s in suffixes) else "fasta"
    return handle, fmt


def write_read_sample(paths, output, maximum):
    from Bio import SeqIO

    paths = [Path(p) for p in paths]
    per_file = max(1, math.ceil(maximum / len(paths)))
    count = 0
    with open(output, "w") as out:
        for file_no, path in enumerate(paths):
            handle, fmt = open_reads(path)
            with handle:
                for read_no, rec in enumerate(SeqIO.parse(handle, fmt)):
                    if read_no >= per_file or count >= maximum:
                        break
                    rec.id = f"r{file_no}_{read_no}"
                    rec.name = rec.id
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
                    count += 1
    return count


def rank_with_reads(read_fasta, subject_fasta, lengths):
    cmd = [
        "blastn", "-task", "blastn", "-query", str(read_fasta),
        "-subject", str(subject_fasta), "-evalue", "1e-8", "-max_hsps", "1",
        "-outfmt", "6 qseqid sseqid sstart send pident length bitscore",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "blastn failed")
    aggregate = {}
    for line in result.stdout.splitlines():
        qid, sid, ss, se, pid, aln_len, bits = line.split("\t")
        row = aggregate.setdefault(
            sid, {"reads": set(), "intervals": [], "aligned_bp": 0,
                  "weighted_pid": 0.0, "bits": 0.0})
        lo, hi = sorted((int(ss), int(se)))
        aln_len = int(aln_len)
        row["reads"].add(qid)
        row["intervals"].append((lo, hi))
        row["aligned_bp"] += aln_len
        row["weighted_pid"] += float(pid) * aln_len
        row["bits"] += float(bits)
    ranking = []
    for sid, row in aggregate.items():
        ref_len = lengths.get(sid) or 1
        coverage = union_coverage(row["intervals"]) / ref_len
        identity = row["weighted_pid"] / row["aligned_bp"] if row["aligned_bp"] else 0.0
        ranking.append((sid, coverage, len(row["reads"]), row["bits"], identity))
    return sorted(ranking, key=lambda x: (x[1], x[2], x[3], x[4]), reverse=True)


def main():
    from Bio import SeqIO

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", nargs="+", required=True, type=Path)
    ap.add_argument("--refdb-dir", required=True, type=Path)
    ap.add_argument("--group", required=True)
    ap.add_argument("--taxon", default="")
    ap.add_argument("--family", default="")
    ap.add_argument("--order", default="")
    ap.add_argument("--class-name", default="")
    ap.add_argument("--top-n", type=int, default=5)
    ap.add_argument("--read-subsample", type=int, default=50000)
    ap.add_argument("--max-candidates", type=int, default=50)
    ap.add_argument("--out-seed-fasta", required=True, type=Path)
    ap.add_argument("--out-label-fasta", required=True, type=Path)
    ap.add_argument("--out-gb", required=True, type=Path)
    ap.add_argument("--out-status", required=True, type=Path)
    args = ap.parse_args()

    def finish(state, message):
        args.out_status.write_text(f"{state}\t{message}\n")
        print(f"[select_fallback_seed] {state}: {message}", file=sys.stderr)
        raise SystemExit(0)

    db_fasta, manifest_path, _features_path = refdb_record.refdb_paths(
        args.refdb_dir, args.group)
    label_path = Path(str(db_fasta).replace(".fasta", ".label.fasta"))
    if not db_fasta.exists() or not manifest_path.exists() or not label_path.exists():
        finish("NONE", f"incomplete reference database for {args.group}")

    rows = refdb_record.read_manifest(manifest_path)
    candidates, tier, missed = taxonomy_shortlist(
        rows, args.taxon, args.family, args.order, args.class_name,
        max(1, args.max_candidates))
    if not candidates:
        finish("NONE", f"no records in {args.group} manifest")

    wanted = set(candidates)
    records = {rec.id: rec for rec in SeqIO.parse(str(db_fasta), "fasta")
               if rec.id in wanted}
    candidates = [acc for acc in candidates if acc in records]
    if not candidates:
        finish("NONE", f"no shortlisted records materialised for {args.group}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        reads_fasta = tmpdir / "reads.fasta"
        subject_fasta = tmpdir / "candidates.fasta"
        n_reads = write_read_sample(args.reads, reads_fasta, max(1, args.read_subsample))
        SeqIO.write([records[a] for a in candidates], str(subject_fasta), "fasta")
        try:
            ranked = rank_with_reads(
                reads_fasta, subject_fasta,
                {a: len(records[a].seq) for a in candidates})
        except Exception as exc:
            finish("NONE", f"read ranking failed: {exc}")

    ranked_ids = [row[0] for row in ranked]
    chosen = (ranked_ids or candidates)[:max(1, args.top_n)]
    SeqIO.write([records[a] for a in chosen], str(args.out_seed_fasta), "fasta")
    n_labels = subset_label_db(label_path, chosen, args.out_label_fasta)
    reference = refdb_record.materialize(args.refdb_dir, args.group, chosen[0])
    if not n_labels or reference is None:
        finish("NONE", f"could not materialise labels/reference for {','.join(chosen)}")
    SeqIO.write([reference], str(args.out_gb), "genbank")

    gap = ",".join(missed) if missed else "none"
    if ranked:
        best = ranked[0]
        support = (f"best={best[0]} ref_cov={best[1]:.3f} "
                   f"reads={best[2]} pid={best[4]:.1f}")
    else:
        support = f"no_read_hits fallback={chosen[0]}"
    finish(
        "SELECTED_SEED",
        f"fallback=taxonomy_reads tier={tier} reference_gap={gap} "
        f"candidates={len(candidates)} sampled_reads={n_reads} "
        f"chosen={','.join(chosen)} {support}")


if __name__ == "__main__":
    main()
