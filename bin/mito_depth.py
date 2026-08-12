#!/usr/bin/env python3
"""Uniform, cross-platform mitogenome read depth.

Every assembler in this pipeline used to report a different quantity under the
name "coverage", so the numbers were never comparable:

  GetOrganelle  k-mer coverage off its assembly graph (~0.2x true depth), and
                measured on the reduced read set its --reduce-reads-for-coverage
                default selects.
  MitoHiFi      per-base depth of the reads it recruited *by mapping to a
                related-species reference*, so a divergent reference silently
                depresses the number.
  Oatk          nothing at all.

This script replaces all three with one definition: **mean per-base depth of the
sample's own reads remapped to the assembly that actually goes to annotation**.
Remapping to self rather than to a reference is the whole point -- a genuinely
divergent mitogenome is measured just as accurately as a well-referenced one,
which is exactly the bias that made the MitoHiFi number untrustworthy.

Three things it gets right that a naive `samtools depth` does not:

1. **Circular molecules are folded.** A closed mitogenome is stored linearised,
   so reads spanning the origin clip and depth collapses at both ends. Mapping to
   a head-to-tail doubled reference and folding (d[p] += d[p+L]) recovers the true
   circular profile. Under --secondary=no each read has exactly one primary
   placement, so folding maps every covered base to exactly one monomer
   coordinate and cannot double count.
2. **NUMTs are filtered by gap-compressed identity, not raw NM/aligned-length.**
   A true read spanning a 30 bp control-region indel has NM=30 on a 150 bp read;
   raw identity calls that 0.75 and discards it, punching a depth hole in the
   D-loop of exactly the VNTR-bearing samples we care about. Gap-compressed
   identity scores it 0.992 and keeps it. minimap2's own de:f: tag is preferred
   when present and the CIGAR+NM formula reproduces it to 4 decimal places.
3. **MAPQ is deliberately ignored.** On a doubled reference every read maps to two
   equally good places, so MAPQ is 0 by construction; and even undoubled there is
   only one reference sequence, so nothing competes. MAPQ carries no information
   here. Identity and aligned-fraction do.

Known limits, stated rather than hidden:
  - Recent NUMTs (<5% diverged) pass the short-read filter. Unavoidable at 150 bp.
    They are normally <1% of mito depth, but that argument is weakest exactly
    where the metric matters most (HiC / low-mito-content libraries at 5-20x), so
    reads_fail_identity is reported rather than silently dropped.
  - Depth over an uncollapsed tandem repeat is diluted by 1/copies, inherent to
    counting each read once.
  - Deletions do not contribute depth, matching samtools depth semantics.

Always exits 0. On any failure it writes a header-only TSV so the downstream
joins in the workflow never break on a missing file.

NOTE: this runs in the MitoHiFi container, whose python3 is 3.6.9. Keep it
3.6-compatible: no subprocess capture_output=/text= (3.7+), no
Path.unlink(missing_ok=) (3.8+), no statistics.fmean (3.8+).
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import zlib
from array import array

# Order matters: this is the header written to <prefix>.mito_depth.tsv.
#
# Column names deliberately avoid the substrings "rejected", "divergent",
# "low_coverage", "low confidence" and "nuclear mitochondrial":
# mitogenome_assembly_summary.py's has_numt_signal() greps the *text* of every
# .tsv in a run for those, so a column called e.g. reads_rejected would flag
# every assembly in the cohort as NUMT-contaminated.
#
# mean_coverage / coverage_cv are intentional duplicates of mean_depth /
# depth_cv so the summary's existing first_numeric(row, ["mean_coverage", ...])
# lookup picks this file up with no special casing.
COLUMNS = [
    "sample",
    "target_fasta",
    "target_length_bp",
    "n_contigs",
    "circular_doubled",
    "preset",
    "sequencing_type",
    "mean_depth",
    "median_depth",
    "sd_depth",
    "depth_cv",
    "p10_depth",
    "min_depth",
    "breadth_1x",
    "breadth_10x",
    "breadth_20x",
    "mito_mapped_reads",
    "reads_fail_identity",
    "reads_fail_clip",
    "supplementary_dropped",
    "total_reads",
    "mito_read_fraction",
    "mean_identity",
    "min_identity",
    "min_aligned_frac",
    "subsampled",
    "subsample_fraction",
    "scale_factor",
    "depth_method",
    "mean_coverage",
    "coverage_cv",
]

DEPTH_METHOD = "remap_full_v1"

# mawk 1.3.3 is the awk in the MitoHiFi container and has no bitwise and(), so
# the SAM flag tests are arithmetic. This prefilter exists purely for speed: it
# discards the ~99.9% of a WGS library that does not touch the mitogenome at C
# speed, so python only ever parses real mito alignments. Survivors are passed
# through as untouched SAM lines, which keeps a single parser in python for both
# this path and the --sam test path.
TRIAGE_AWK = r"""
/^@/ { next }
{
    f = $2
    sec = int(f / 256) % 2
    sup = int(f / 2048) % 2
    unm = int(f / 4) % 2
    if (sec == 0 && sup == 0) total++
    if (unm == 1 || sec == 1 || sup == 1) {
        if (sup == 1) supp++
        next
    }
    print
}
END { printf "#TOTAL\t%d\t%d\n", total, supp }
"""

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# CIGAR ops that consume the reference and represent an aligned base of the read.
REF_COVERING = ("M", "=", "X")
# Consume the reference but contribute no depth (samtools depth semantics).
REF_SKIPPING = ("D", "N")
# Consume the query only.
QUERY_ONLY = ("I", "S")


def open_maybe_gzip(path):
    """Same helper as bin/rank_reference_candidates.py; 3.6-safe."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def remove_quietly(path):
    """Path.unlink(missing_ok=True) is 3.8+; this container is on 3.6."""
    try:
        os.remove(str(path))
    except OSError:
        pass


def read_fasta(path):
    """Return [(name, seq), ...]. Raises on an unreadable/empty file."""
    records = []
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise ValueError("no FASTA records in {}".format(path))
    return records


def write_doubled(records, out_path):
    """Head-to-tail doubled reference, mirroring bin/check_circularity.py."""
    with open(out_path, "w") as handle:
        for name, seq in records:
            handle.write(">{}_doubled\n".format(name))
            doubled = seq + seq
            for i in range(0, len(doubled), 80):
                handle.write(doubled[i:i + 80] + "\n")


def should_double(circular, n_records, fasta_stem):
    """Double only a molecule we actually believe is one closed circle.

    The _concat guard is load-bearing: this script runs on the SANITISE_FASTA
    output, which collapses every multi-contig assembly into a single record, so
    n_records == 1 no longer distinguishes a real circle from a blind
    concatenation of scaffolds. Doubling a _concat molecule would fabricate a
    junction that does not exist.
    """
    if circular is not True:
        return False
    if n_records != 1:
        return False
    if fasta_stem.endswith("_concat"):
        return False
    return True


def parse_circular(value):
    if value is None:
        return None
    token = str(value).strip().lower()
    if token in ("true", "t", "yes", "y", "1"):
        return True
    if token in ("false", "f", "no", "n", "0"):
        return False
    return None


def cigar_blocks(cigar):
    """Yield (op, length) pairs; '*' yields nothing."""
    if not cigar or cigar == "*":
        return
    for length, op in CIGAR_RE.findall(cigar):
        yield op, int(length)


def alignment_geometry(cigar):
    """Return (ref_blocks, ref_span, aligned_query, query_len, indel_events,
    indel_bases) for one CIGAR.

    ref_blocks are (start_offset, end_offset) pairs relative to the alignment
    start, covering only the reference bases the read actually places a base on.
    """
    ref_blocks = []
    ref_pos = 0
    aligned_query = 0        # M/=/X only: the identity denominator
    query_aligned_span = 0   # M/=/X/I: query bases participating in the alignment
    query_len = 0            # everything, including clips: the full read length
    indel_events = 0
    indel_bases = 0
    leading_clip = 0
    trailing_clip = 0
    seen_aligned = False
    for op, length in cigar_blocks(cigar):
        if op in REF_COVERING:
            ref_blocks.append((ref_pos, ref_pos + length))
            ref_pos += length
            aligned_query += length
            query_aligned_span += length
            query_len += length
            seen_aligned = True
            trailing_clip = 0
        elif op in REF_SKIPPING:
            if op == "D":
                indel_events += 1
                indel_bases += length
            ref_pos += length
        elif op in QUERY_ONLY:
            if op == "I":
                indel_events += 1
                indel_bases += length
                # An inserted base IS part of the alignment, just not of the
                # reference. Excluding it here would make a fully-aligned read
                # carrying a large insertion (e.g. 50M100I50M) look 50% clipped and
                # fail the aligned-fraction floor -- exactly the large-indel case
                # the gap-compressed identity is designed to keep.
                query_aligned_span += length
                query_len += length
            else:  # S, soft clip: part of the read, not aligned
                query_len += length
                if not seen_aligned:
                    leading_clip += length
                else:
                    trailing_clip += length
        elif op == "H":
            # Hard clip: not present in SEQ but still part of the original read,
            # so it belongs in the aligned-fraction denominator.
            query_len += length
            if not seen_aligned:
                leading_clip += length
            else:
                trailing_clip += length
    return {
        "ref_blocks": ref_blocks,
        "ref_span": ref_pos,
        "aligned_query": aligned_query,
        "query_aligned_span": query_aligned_span,
        "query_len": query_len,
        "indel_events": indel_events,
        "indel_bases": indel_bases,
        "leading_clip": leading_clip,
        "trailing_clip": trailing_clip,
    }


def gap_compressed_identity(nm, aligned_query, indel_events, indel_bases):
    """Identity counting each indel as ONE difference regardless of its length.

    minimap2's NM is mismatches + inserted + deleted bases, so the mismatch count
    is NM minus the gap bases. Verified against minimap2 2.24's own de:f: tag:
    a 150 bp read with 74M30D46M / NM:i:30 gives 1 - 1/121 = 0.9917, and
    minimap2 reports de:f:0.0083.
    """
    denominator = aligned_query + indel_events
    if denominator <= 0:
        return None
    mismatches = nm - indel_bases
    if mismatches < 0:
        mismatches = 0
    return 1.0 - float(mismatches + indel_events) / denominator


def parse_tags(fields):
    """Pull NM:i: and de:f: out of the optional SAM columns."""
    nm = None
    de = None
    for field in fields:
        if field.startswith("NM:i:"):
            try:
                nm = int(field[5:])
            except ValueError:
                nm = None
        elif field.startswith("de:f:"):
            try:
                de = float(field[5:])
            except ValueError:
                de = None
    return nm, de


class DepthAccumulator(object):
    """Difference array, so each aligned block is O(1) rather than O(length).

    Prefix-summed once at the end. A 1000x mitogenome would otherwise be ~17M
    per-base increments; this makes it a few hundred thousand.
    """

    def __init__(self, length):
        self.length = length
        self.diff = array("l", [0] * (length + 1))

    def add(self, start, end):
        if end <= 0 or start >= self.length:
            return
        if start < 0:
            start = 0
        if end > self.length:
            end = self.length
        self.diff[start] += 1
        self.diff[end] -= 1

    def depths(self):
        out = array("l", [0] * self.length)
        running = 0
        for i in range(self.length):
            running += self.diff[i]
            out[i] = running
        return out


def fold_circular(depths, monomer_length):
    """d[p] += d[p + L]; exact because --secondary=no places each read once."""
    folded = array("l", [0] * monomer_length)
    for i in range(monomer_length):
        folded[i] = depths[i]
    for i in range(monomer_length, len(depths)):
        folded[i - monomer_length] += depths[i]
    return folded


def summarise(depths):
    """mean / median / sd / cv / p10 / min / breadth, all in one pass-ish."""
    n = len(depths)
    if n == 0:
        return {}
    total = 0
    for value in depths:
        total += value
    mean = float(total) / n

    variance_total = 0.0
    for value in depths:
        delta = value - mean
        variance_total += delta * delta
    sd = (variance_total / n) ** 0.5

    ordered = sorted(depths)
    if n % 2:
        median = float(ordered[n // 2])
    else:
        median = (ordered[n // 2 - 1] + ordered[n // 2]) / 2.0
    p10 = float(ordered[max(0, int(0.10 * n) - 1)]) if n else 0.0

    at_least = {}
    for threshold in (1, 10, 20):
        count = 0
        for value in depths:
            if value >= threshold:
                count += 1
        at_least[threshold] = float(count) / n

    return {
        "mean": mean,
        "median": median,
        "sd": sd,
        "cv": (sd / mean) if mean > 0 else None,
        "p10": p10,
        "min": float(ordered[0]),
        "breadth_1x": at_least[1],
        "breadth_10x": at_least[10],
        "breadth_20x": at_least[20],
    }


def subsample_fastq(read_paths, out_path, fraction):
    """Uniform Bernoulli draw over the WHOLE file, keyed on a stable hash.

    Deliberately not the head-of-file draw that rank_reference_candidates.py
    uses: that one is fine for comparative candidate ranking (the bias cancels
    across candidates) but would bias an absolute depth estimate, since reads are
    ordered by flowcell position. crc32 of the read name is stable, uniform, and
    keeps R1/R2 together because a pair shares its name.

    Returns the number of records kept.
    """
    cutoff = int(fraction * (1 << 32))
    kept = 0
    with open(out_path, "w") as out:
        for read_path in read_paths:
            with open_maybe_gzip(read_path) as handle:
                while True:
                    header = handle.readline()
                    if not header:
                        break
                    seq = handle.readline()
                    plus = handle.readline()
                    qual = handle.readline()
                    if not qual:
                        break
                    name = header[1:].split()[0] if len(header) > 1 else ""
                    # Strip a trailing /1 or /2 so mates hash identically.
                    if name.endswith("/1") or name.endswith("/2"):
                        name = name[:-2]
                    if zlib.crc32(name.encode("utf-8")) & 0xFFFFFFFF < cutoff:
                        out.write(header)
                        out.write(seq)
                        out.write(plus)
                        out.write(qual)
                        kept += 1
    return kept


def sam_line_source(args, workdir, reference):
    """Yield SAM lines, either from minimap2 (piped through the awk triage) or
    from a --sam file for testing. Returns (iterator, process-or-None)."""
    if args.sam:
        if args.sam == "-":
            return sys.stdin, None
        return open(args.sam), None

    awk_path = os.path.join(workdir, "triage.awk")
    with open(awk_path, "w") as handle:
        handle.write(TRIAGE_AWK)

    reads = " ".join('"{}"'.format(path) for path in args.reads)
    stderr_path = os.path.join(workdir, "minimap2.stderr")
    # pipefail is essential: without it the pipeline's status is awk's, and awk
    # succeeds happily on empty input. A missing/crashed/OOM-killed minimap2 would
    # then look like "mapped 0 reads" and be written to SQL as a genuine depth of
    # zero. Needs bash explicitly, since /bin/sh may be dash.
    command = (
        'set -o pipefail; minimap2 -ax {preset} --secondary=no -t {threads} '
        '"{ref}" {reads} 2> "{err}" | awk -f "{awk}"'
    ).format(
        preset=args.preset,
        threads=args.threads,
        ref=reference,
        reads=reads,
        err=stderr_path,
        awk=awk_path,
    )
    process = subprocess.Popen(
        command, shell=True, executable="/bin/bash",
        stdout=subprocess.PIPE, universal_newlines=True,
    )
    process.mito_depth_stderr = stderr_path
    return process.stdout, process


def write_tsv(out_path, row):
    with open(out_path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        if row is not None:
            handle.write("\t".join(str(row.get(column, "")) for column in COLUMNS) + "\n")


def fmt(value, places=6):
    if value is None:
        return "NA"
    return "{:.{}f}".format(value, places).rstrip("0").rstrip(".") or "0"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Assembly that goes to annotation.")
    parser.add_argument("--reads", nargs="*", default=[], help="FASTQ(.gz) read file(s).")
    parser.add_argument("--sample", required=True, help="meta.mt_assembly_prefix (DB identity).")
    parser.add_argument("--out", required=True, help="Output <prefix>.mito_depth.tsv.")
    parser.add_argument("--preset", default="map-hifi", choices=["sr", "map-hifi"])
    parser.add_argument("--sequencing-type", default="")
    parser.add_argument("--circular", default=None, help="true / false / null.")
    parser.add_argument("--min-identity", type=float, default=0.95)
    parser.add_argument("--min-aligned-frac", type=float, default=0.80)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--subsample-fraction",
        type=float,
        default=0.0,
        help="0 = use every read. 0 < f <= 1 keeps a uniform random fraction.",
    )
    parser.add_argument("--sam", default=None, help="Read SAM from here instead of mapping (tests).")
    parser.add_argument("--workdir", default=".")
    args = parser.parse_args()

    # The label must name the molecule that is about to be measured. This is checked BEFORE
    # the fail-open wrapper below and exits non-zero on purpose: a mismatch here is not a
    # measurement problem to be degraded gracefully, it is the depth being attributed to the
    # wrong assembly, and writing a placeholder row would hide it.
    #
    # This is the exact defect that motivated the check. A reseeded sample was remapped
    # correctly against OG5.hic.250522.getorg1770reseed.fasta, but --sample was handed the
    # sample-level prefix, so the TSV said "OG5.hic.250522.getorg1770" and the resulting
    # 257.8x was filed against the superseded first-pass assembly while the reseed -- the
    # molecule actually annotated and submitted -- recorded no depth at all.
    fasta_stem = re.sub(r"\.(fa|fasta|fna)$", "", os.path.basename(args.fasta))
    if args.sample != fasta_stem:
        print(
            "[mito_depth] --sample '{}' does not name the FASTA being measured ('{}'). "
            "The assembly's identity (meta.mt_assembly_prefix) must equal its FASTA basename; "
            "whichever stage produced this FASTA has not stamped it.".format(
                args.sample, fasta_stem
            ),
            file=sys.stderr,
        )
        return 1

    try:
        run(args)
    except Exception as exc:  # noqa: BLE001 - fail open, never break the workflow
        print("[mito_depth] failed: {}".format(exc), file=sys.stderr)
        try:
            write_tsv(args.out, None)
        except Exception:  # noqa: BLE001
            pass
    return 0


def run(args):
    records = read_fasta(args.fasta)
    monomer_length = sum(len(seq) for _name, seq in records)
    n_records = len(records)
    fasta_stem = re.sub(r"\.(fa|fasta|fna)$", "", os.path.basename(args.fasta))

    circular = parse_circular(args.circular)
    doubled = should_double(circular, n_records, fasta_stem)

    reference = args.fasta
    workdir = args.workdir
    if doubled:
        reference = os.path.join(workdir, "doubled.fasta")
        write_doubled(records, reference)

    reference_length = monomer_length * 2 if doubled else monomer_length

    # Optional unbiased subsample. Off by default: it saves the minimap2 seeding
    # cost but still decompresses every read, so it is not a free win.
    subsample_path = None
    reads = list(args.reads)
    subsampled = False
    if args.subsample_fraction and 0 < args.subsample_fraction < 1 and reads:
        subsample_path = os.path.join(workdir, "subsample.fastq")
        subsample_fastq(reads, subsample_path, args.subsample_fraction)
        reads = [subsample_path]
        subsampled = True
    args.reads = reads

    accumulator = DepthAccumulator(reference_length)
    mapped = 0
    fail_identity = 0
    fail_clip = 0
    supplementary = 0
    total_reads = 0
    identity_total = 0.0
    # With the awk prefilter, unmapped records never reach python, so the total
    # can ONLY come from awk's #TOTAL sentinel. On the --sam test path there is no
    # prefilter and python counts for itself.
    use_awk = args.sam is None
    # Set by the awk END sentinel. Its absence means the stream was truncated, so
    # the run is a failure rather than a legitimate zero.
    saw_total = not use_awk

    source, process = sam_line_source(args, workdir, reference)
    try:
        for line in source:
            if not line:
                continue
            if line.startswith("#TOTAL"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    total_reads = int(parts[1])
                if len(parts) >= 3:
                    supplementary = int(parts[2])
                saw_total = True
                continue
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            secondary = bool(flag & 0x100)
            supplementary_flag = bool(flag & 0x800)
            unmapped = bool(flag & 0x4)
            if not use_awk:
                if not secondary and not supplementary_flag:
                    total_reads += 1
                if supplementary_flag:
                    supplementary += 1
            if unmapped or secondary or supplementary_flag:
                continue

            pos = int(fields[3]) - 1
            cigar = fields[5]
            geometry = alignment_geometry(cigar)
            ref_blocks = geometry["ref_blocks"]
            if not ref_blocks:
                continue

            nm, de = parse_tags(fields[11:])
            if de is not None:
                identity = 1.0 - de
            elif nm is not None:
                identity = gap_compressed_identity(
                    nm, geometry["aligned_query"],
                    geometry["indel_events"], geometry["indel_bases"],
                )
            else:
                identity = None
            if identity is not None and identity < args.min_identity:
                fail_identity += 1
                continue

            # Aligned-fraction floor, with a DIRECTIONAL waiver for clipping the
            # reference itself explains. Only the clip on the side that actually
            # runs off the reference is excused: an alignment starting at
            # position 0 excuses leading clip, one ending at the last base excuses
            # trailing clip. A blanket "at either edge -> skip the test" waiver
            # would admit arbitrary junk, e.g. 1M149S at position 1 (0.7% of the
            # read aligned) would count as a mitochondrial read.
            query_len = geometry["query_len"]
            excused = 0
            if pos <= 0:
                excused += geometry["leading_clip"]
            if (pos + geometry["ref_span"]) >= reference_length:
                excused += geometry["trailing_clip"]
            testable_len = query_len - excused
            if testable_len > 0:
                if float(geometry["query_aligned_span"]) / testable_len < args.min_aligned_frac:
                    fail_clip += 1
                    continue

            for start, end in ref_blocks:
                accumulator.add(pos + start, pos + end)
            mapped += 1
            if identity is not None:
                identity_total += identity
    finally:
        # Only clean up here. The failure checks live below, outside the finally:
        # raising from a finally would replace any exception already propagating
        # out of the loop with a less informative one.
        if hasattr(source, "close") and source is not sys.stdin:
            source.close()
        if process is not None:
            returncode = process.wait()

    if process is not None:
        # A broken mapper must never look like a real measurement. Without these
        # checks a missing, crashed or OOM-killed minimap2 leaves awk emitting
        # "#TOTAL 0 0" and exiting cleanly, which would be written to the database
        # as a genuine mean_depth of 0. Raising lets main() fail open with a
        # header-only TSV instead.
        if returncode != 0:
            detail = ""
            try:
                with open(process.mito_depth_stderr) as handle:
                    detail = " | ".join(handle.read().strip().splitlines()[-5:])
            except (OSError, IndexError):
                pass
            raise RuntimeError(
                "minimap2/awk pipeline exited {}: {}".format(returncode, detail)
            )
        if not saw_total:
            raise RuntimeError(
                "mapping produced no #TOTAL record; treating as a failed run "
                "rather than a depth of zero"
            )

    depths = accumulator.depths()
    if doubled:
        depths = fold_circular(depths, monomer_length)

    # Scale the per-base depths BEFORE summarising, not the summary statistics
    # afterwards. Scaling afterwards would leave breadth_10x / breadth_20x counted
    # against the sampled depths, so a true 20x position sampled at 25% would show
    # as 5x and fail its own breadth threshold.
    scale_factor = 1.0
    realised_fraction = None
    if subsampled and args.subsample_fraction > 0:
        realised_fraction = args.subsample_fraction
        scale_factor = 1.0 / args.subsample_fraction
        depths = [value * scale_factor for value in depths]

    stats = summarise(depths)

    mito_fraction = None
    if total_reads:
        mito_fraction = float(mapped) / total_reads

    row = {
        "sample": args.sample,
        "target_fasta": os.path.basename(args.fasta),
        "target_length_bp": monomer_length,
        "n_contigs": n_records,
        "circular_doubled": "true" if doubled else "false",
        "preset": args.preset,
        "sequencing_type": args.sequencing_type,
        "mean_depth": fmt(stats.get("mean")),
        "median_depth": fmt(stats.get("median")),
        "sd_depth": fmt(stats.get("sd")),
        "depth_cv": fmt(stats.get("cv")),
        "p10_depth": fmt(stats.get("p10")),
        "min_depth": fmt(stats.get("min")),
        "breadth_1x": fmt(stats.get("breadth_1x")),
        "breadth_10x": fmt(stats.get("breadth_10x")),
        "breadth_20x": fmt(stats.get("breadth_20x")),
        "mito_mapped_reads": mapped,
        "reads_fail_identity": fail_identity,
        "reads_fail_clip": fail_clip,
        "supplementary_dropped": supplementary,
        "total_reads": total_reads,
        "mito_read_fraction": fmt(mito_fraction, 8),
        "mean_identity": fmt(identity_total / mapped) if mapped else "NA",
        "min_identity": fmt(args.min_identity),
        "min_aligned_frac": fmt(args.min_aligned_frac),
        "subsampled": "true" if subsampled else "false",
        "subsample_fraction": fmt(realised_fraction) if realised_fraction else "NA",
        "scale_factor": fmt(scale_factor),
        "depth_method": DEPTH_METHOD,
        "mean_coverage": fmt(stats.get("mean")),
        "coverage_cv": fmt(stats.get("cv")),
    }
    write_tsv(args.out, row)

    if subsample_path:
        remove_quietly(subsample_path)
    if doubled and reference != args.fasta:
        remove_quietly(reference)

    print(
        "[mito_depth] {} mean_depth={} cv={} breadth_1x={} mapped={}/{}".format(
            args.sample,
            row["mean_depth"],
            row["depth_cv"],
            row["breadth_1x"],
            mapped,
            total_reads,
        )
    )


if __name__ == "__main__":
    sys.exit(main())
