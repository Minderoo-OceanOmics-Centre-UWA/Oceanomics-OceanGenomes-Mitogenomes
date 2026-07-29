#!/usr/bin/env python3
"""Choose the best mitogenome reference for a sample from several candidates, by
mapping a subsample of the sample's own reads against each one.

Why read evidence rather than taxonomy: MITOHIFI_FINDMITOREFERENCE returns the
first complete mitogenome it meets walking up the NCBI lineage, so for a taxon with
no congeneric record the reference is an arbitrary member of whatever rank the walk
reached. MitoHiFi then builds the assembly by recruiting reads that map to that
reference, so "how many of this sample's reads does the candidate actually recruit"
is precisely the quantity that determines whether the assembly succeeds -- and it
is measurable before assembly, from reads alone.

Taxonomy cannot answer this on its own: NCBI lineages are uneven, and two
same-rank candidates can differ enormously in how well they recruit. Sequence
evidence from the sample settles it directly.

Scoring: each candidate is scored by the number of read bases that align to it
(PAF residue matches), normalised by candidate length so a longer reference does
not win merely by being longer. The candidate with the highest score wins; ties
keep the earlier candidate, which is the one findMitoReference ranked closest
taxonomically, so this can only improve on the current behaviour.

Writes <prefix>.reference_ranking.tsv (every candidate and its score, so the choice
is auditable) and copies the winner to <prefix>.chosen_reference.{fasta,gb}.
Always exits 0: on any failure it falls back to the first candidate, which is what
the pipeline would have used anyway.

Usage:
    rank_reference_candidates.py --candidate-dir candidates --reads r1.fq.gz [r2.fq.gz] \
        --preset map-hifi --prefix OG123.hifi.240101.v323mitohifi [--subsample-reads 50000]
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

# NOTE: this script runs in the MitoHiFi container, whose python3 is 3.6.9. Keep it
# 3.6-compatible: no subprocess capture_output=/text= (3.7+), no Path.unlink(
# missing_ok=) (3.8+). bin/check_getorganelle.py uses the same 3.6-safe idiom.


def remove_quietly(path):
    """Path.unlink(missing_ok=True) is 3.8+; this container is on 3.6."""
    try:
        os.remove(str(path))
    except OSError:
        pass


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def subsample_reads(read_paths, out_fastq, n_reads):
    """Write the first n_reads records across the given FASTQ(s) to out_fastq.

    Head-of-file rather than a random draw: it needs no index and no second pass
    over what can be tens of GB of HiFi reads, and every candidate is scored on the
    same reads, so any positional bias applies equally to all of them and cancels
    out of the comparison.
    """
    written = 0
    with open(out_fastq, "w") as out:
        for path in read_paths:
            if written >= n_reads:
                break
            try:
                with open_maybe_gzip(path) as fh:
                    for line_no, line in enumerate(fh):
                        out.write(line)
                        if line_no % 4 == 3:
                            written += 1
                            if written >= n_reads:
                                break
            except OSError as exc:
                print(f"[rank_reference] could not read {path}: {exc}", file=sys.stderr)
    return written


def score_candidate(candidate_fasta, fastq, preset, threads):
    """Aligned read bases per reference base, from minimap2 PAF residue matches."""
    try:
        proc = subprocess.run(
            ["minimap2", "-x", preset, "-t", str(threads), "--secondary=no",
             str(candidate_fasta), str(fastq)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except (subprocess.CalledProcessError, OSError) as exc:
        print(f"[rank_reference] minimap2 failed on {candidate_fasta.name}: {exc}",
              file=sys.stderr)
        return 0.0, 0, 0

    matches, reads = 0, set()
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        reads.add(fields[0])
        matches += int(fields[9])       # PAF col 10: residue matches

    ref_len = candidate_length(candidate_fasta)
    return (matches / ref_len if ref_len else 0.0), matches, len(reads)


def candidate_length(fasta):
    total = 0
    try:
        with open(fasta) as fh:
            for line in fh:
                if not line.startswith(">"):
                    total += len(line.strip())
    except OSError:
        return 0
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate-dir", required=True, type=Path)
    ap.add_argument("--reads", required=True, nargs="+", type=Path)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--preset", default="map-hifi",
                    help="minimap2 preset: map-hifi for HiFi, sr for short reads.")
    ap.add_argument("--subsample-reads", type=int, default=50000)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    ranking_out = Path(f"{args.prefix}.reference_ranking.tsv")
    chosen_fa = Path(f"{args.prefix}.chosen_reference.fasta")
    chosen_gb = Path(f"{args.prefix}.chosen_reference.gb")

    # A candidate is usable only as a fasta+gb pair: MitoHiFi needs both.
    candidates = []
    for fasta in sorted(args.candidate_dir.glob("*.fasta")):
        gb = fasta.with_suffix(".gb")
        if gb.exists() and fasta.stat().st_size > 0 and gb.stat().st_size > 0:
            candidates.append((fasta, gb))

    if not candidates:
        ranking_out.write_text("accession\tscore\tmatched_bases\tmapped_reads\tchosen\n")
        print("[rank_reference] no usable candidates; leaving reference unchanged",
              file=sys.stderr)
        sys.exit(0)

    def emit(winner_idx, rows):
        fasta, gb = candidates[winner_idx]
        shutil.copyfile(fasta, chosen_fa)
        shutil.copyfile(gb, chosen_gb)
        with open(ranking_out, "w") as out:
            out.write("accession\tscore\tmatched_bases\tmapped_reads\tchosen\n")
            for i, (acc, score, matched, mapped) in enumerate(rows):
                out.write(f"{acc}\t{score:.4f}\t{matched}\t{mapped}\t"
                          f"{'yes' if i == winner_idx else 'no'}\n")
        print(f"[rank_reference] chose {fasta.stem} of {len(candidates)} candidates",
              file=sys.stderr)
        sys.exit(0)

    # Single candidate: nothing to choose between, skip the mapping entirely.
    if len(candidates) == 1:
        emit(0, [(candidates[0][0].stem, 0.0, 0, 0)])

    subsample = Path("reference_rank_subsample.fastq")
    n_reads = subsample_reads(args.reads, subsample, args.subsample_reads)
    if n_reads == 0:
        print("[rank_reference] no reads subsampled; keeping first candidate",
              file=sys.stderr)
        emit(0, [(fa.stem, 0.0, 0, 0) for fa, _ in candidates])

    rows, best_idx, best_score = [], 0, -1.0
    for i, (fasta, _gb) in enumerate(candidates):
        score, matched, mapped = score_candidate(fasta, subsample, args.preset, args.threads)
        rows.append((fasta.stem, score, matched, mapped))
        # Strict > keeps the earliest (taxonomically closest) candidate on a tie.
        if score > best_score:
            best_idx, best_score = i, score

    remove_quietly(subsample)
    emit(best_idx, rows)


if __name__ == "__main__":
    main()
