#!/usr/bin/env python3
"""Re-test circularity for a MitoHiFi assembly and correct was_circular in place.

MitoHiFi's own circularisation check (terminal-overlap detection on the trimmed
final contig) frequently yields false negatives on hifiasm assemblies that are
genuinely circular: the assembler closes the unitig but, after MitoHiFi rotates
and trims it, there is no longer a self-overlap for MitoHiFi to re-detect, so it
writes ``was_circular = False``. That false flag then propagates to the SQL
upload (recorded as a "scaffold" instead of "circular genome") and to the
assembly summary (sample flagged for manual review with reason
``not_circularised``).

This script applies two independent circularity signals to a non-circular
MitoHiFi result:

  1. Read-spanning junction test (primary, assembler-agnostic): the HiFi reads
     MitoHiFi already mapped to the final mitogenome are remapped to a *doubled*
     reference (sequence concatenated to itself). A genuinely circular molecule
     produces reads that bridge the artificial linearisation point at position
     ``L``; a linear molecule does not. We count reads whose alignment extends at
     least ``--min-overhang`` bp either side of the junction.

  2. hifiasm c/l contig flag (corroboration): hifiasm names circular unitigs with
     a trailing ``c`` and linear ones with ``l``. The contig that became the final
     mitogenome is read from the contig-stats table; a ``c`` suffix is independent
     evidence of closure.

If either signal indicates circularity, ``was_circular`` is rewritten to ``True``
in a copy of the contig-stats table (the original MitoHiFi verdict and all
supporting evidence are preserved in a separate ``*.circularity_check.tsv``). The
script always writes its outputs and exits 0, even when the remap tooling fails,
so it can never drop a sample from the downstream channels; in that case it falls
back to the hifiasm flag alone and records the failure in the evidence file.

The script also runs a length / tandem-repeat assessment on every assembly,
independent of the circularity verdict. Some mitogenomes (notably many fish)
carry a tandem-repeat array (VNTR) or a duplicated block in the control region
(D-loop); hifiasm assembles one representation of it, often the longest, which
inflates the assembly well beyond the reference and shows up as a repeat region
rather than an assembly overlap. We compare the assembly length to the related
reference (parsed from the contig-stats header), self-align the assembly to
locate any tandem-repeat array, and use the final GenBank gene coordinates to
decide whether that array falls in the control region. When an over-length
assembly coincides with a control-region tandem repeat, the evidence file flags
the sample for manual curation and states the redundant span suggested for
collapse. This is a flag for review, never an automatic edit: the repeat is real
sequence, and the "correct" copy number is a curation decision (length
heteroplasmy means there is no single true length).
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile


def log(msg):
    print(f"[check_circularity] {msg}", file=sys.stderr)


def read_fasta(path):
    """Return (header, sequence) for the first (and expected only) record."""
    header = None
    seq_parts = []
    n_records = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                n_records += 1
                if n_records > 1:
                    break
                header = line[1:].split()[0]
            elif header is not None:
                seq_parts.append(line.strip())
    if n_records > 1:
        log(f"WARNING: {path} has >1 record; using the first ({header}).")
    return header, "".join(seq_parts)


# CIGAR operations that consume the reference.
_REF_CONSUMING = set("MDN=X")


def ref_aligned_length(cigar):
    """Reference span of an alignment from its CIGAR string."""
    total = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in _REF_CONSUMING and num:
                total += int(num)
            num = ""
    return total


def run(cmd, **kwargs):
    log("$ " + (cmd if isinstance(cmd, str) else " ".join(cmd)))
    return subprocess.run(cmd, check=True, **kwargs)


def parse_reference_length(stats_path):
    """Pull the related-mitogenome length from the contig-stats header comment,
    e.g. '# Related mitogenome is 17267 bp long and has 37 genes'."""
    try:
        with open(stats_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    m = re.search(r"is\s+(\d+)\s*bp", line)
                    if m:
                        return int(m.group(1))
                else:
                    break
    except Exception:
        pass
    return None


def _coords_from_location(loc):
    """All integers in a GenBank location string (handles complement/join)."""
    return [int(n) for n in re.findall(r"\d+", loc)]


def parse_gene_intervals(gb_path):
    """Return merged [start, end] intervals (1-based) covered by coding/RNA
    features in a GenBank file, parsed without Biopython."""
    intervals = []
    feat_re = re.compile(r"^\s{3,}(gene|CDS|tRNA|rRNA|mRNA)\s+(\S.*)$")
    try:
        with open(gb_path) as fh:
            for line in fh:
                m = feat_re.match(line.rstrip("\n"))
                if not m:
                    continue
                coords = _coords_from_location(m.group(2))
                if coords:
                    intervals.append((min(coords), max(coords)))
    except Exception as e:
        log(f"Could not parse GenBank {gb_path}: {e}")
        return []
    intervals.sort()
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def control_region(gene_intervals, L):
    """Largest non-coding stretch on the circular molecule -> (start, end).
    Returns None if it cannot be determined."""
    if not gene_intervals or L <= 0:
        return None
    gaps = []
    # Internal gaps between consecutive merged gene intervals.
    for (s1, e1), (s2, e2) in zip(gene_intervals, gene_intervals[1:]):
        if s2 - e1 > 1:
            gaps.append((e1 + 1, s2 - 1))
    # Wrap-around gap (end of last gene -> origin -> start of first gene).
    first_start = gene_intervals[0][0]
    last_end = gene_intervals[-1][1]
    wrap_len = (L - last_end) + (first_start - 1)
    if wrap_len > 0:
        # Represent the wrap as the larger of its two physical arms so the
        # interval we report is contiguous and easy to act on.
        tail = (last_end + 1, L) if last_end < L else None
        head = (1, first_start - 1) if first_start > 1 else None
        for arm in (tail, head):
            if arm:
                gaps.append(arm)
    if not gaps:
        return None
    return max(gaps, key=lambda g: g[1] - g[0])


def find_tandem_repeat(fasta_path, workdir, min_len, min_pident):
    """Self-align the assembly to locate a forward tandem-repeat array.

    Returns a dict (array_start, array_end, array_len, unit_bp, copies, identity)
    or None if no qualifying array is found. The fundamental repeat unit is the
    smallest positive query/subject offset among off-diagonal self-hits; a tandem
    array of unit U self-aligns at offsets U, 2U, 3U..."""
    db = os.path.join(workdir, "selfdb")
    run(
        ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl", "-out", db],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    res = subprocess.run(
        ["blastn", "-query", fasta_path, "-db", db, "-dust", "no",
         "-outfmt", "6 qstart qend sstart send length pident"],
        check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    hits = []
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        qs, qe, ss, se = (int(parts[i]) for i in range(4))
        ln = int(parts[4])
        pid = float(parts[5])
        if ln < min_len or pid < min_pident:
            continue
        # Forward repeats only: subject copy downstream of the query copy.
        if se <= ss or ss <= qs:
            continue
        hits.append((qs, qe, ss, se, ln, pid))
    if not hits:
        return None
    start = min(min(h[0], h[2]) for h in hits)
    end = max(max(h[1], h[3]) for h in hits)
    unit = min(h[2] - h[0] for h in hits)  # smallest offset = fundamental period
    array_len = end - start + 1
    # Require a genuine array (at least ~2 units long) to avoid calling a single
    # short incidental self-match a tandem repeat.
    if unit <= 0 or array_len < max(2 * unit, min_len):
        return None
    total_w = sum(h[4] for h in hits)
    identity = sum(h[4] * h[5] for h in hits) / total_w if total_w else 0.0
    return {
        "array_start": start,
        "array_end": end,
        "array_len": array_len,
        "unit_bp": unit,
        "copies": round(array_len / unit, 1),
        "identity": round(identity, 1),
    }


def assess_length_and_repeats(seq, ref_len, gb_path, workdir, args):
    """Length-vs-reference + tandem-repeat-in-control-region assessment.
    Returns a dict of evidence fields plus a human-readable curation suggestion."""
    L = len(seq)
    out = {
        "assembly_length": L,
        "reference_length": ref_len if ref_len else "NA",
        "length_ratio": "NA",
        "excess_bp": "NA",
        "length_anomaly": "NA",
        "tandem_repeat": "no",
        "repeat_region": "NA",
        "repeat_unit_bp": "NA",
        "repeat_copies": "NA",
        "repeat_identity": "NA",
        "repeat_in_control_region": "NA",
        "anomaly_type": "none",
        "suggested_trim_region": "NA",
        "curation_suggestion": "none",
    }

    length_anomaly = None
    excess = None
    if ref_len:
        ratio = L / ref_len
        excess = L - ref_len
        length_anomaly = ratio >= args.max_length_ratio
        out["length_ratio"] = round(ratio, 3)
        out["excess_bp"] = excess
        out["length_anomaly"] = "yes" if length_anomaly else "no"

    repeat = None
    try:
        repeat = find_tandem_repeat(seq_fasta(seq, workdir), workdir,
                                    args.min_repeat_length, args.min_repeat_identity)
    except subprocess.CalledProcessError as e:
        log(f"Self-alignment failed ({e}); skipping repeat assessment.")
    except Exception as e:  # pragma: no cover - defensive
        log(f"Repeat assessment error ({e}); skipping.")

    cr = None
    if gb_path and os.path.exists(gb_path) and os.path.getsize(gb_path) > 0:
        cr = control_region(parse_gene_intervals(gb_path), L)

    if repeat:
        out["tandem_repeat"] = "yes"
        out["repeat_region"] = f"{repeat['array_start']}-{repeat['array_end']}"
        out["repeat_unit_bp"] = repeat["unit_bp"]
        out["repeat_copies"] = repeat["copies"]
        out["repeat_identity"] = repeat["identity"]
        if cr:
            overlap = min(repeat["array_end"], cr[1]) - max(repeat["array_start"], cr[0])
            out["repeat_in_control_region"] = "yes" if overlap > 0 else "no"

    # Classify the anomaly so the curation action is unambiguous:
    #   - concatemer: the repeat array spans most of the assembly and the length
    #     is a near-integer multiple of the reference -> a tandem genome
    #     duplication (head-to-tail dimer/multimer). Collapse to a monomer.
    #   - control_region_repeat: a localised tandem array inside the D-loop -> a
    #     VNTR / duplicated control region. Real sequence; collapse to a
    #     representative copy or keep and annotate (length heteroplasmy).
    #   - unresolved: over-length with no clean tandem repeat -> manual review for
    #     a partial duplication or NUMT.
    in_cr = out["repeat_in_control_region"]
    array_frac = repeat["array_len"] / L if repeat else 0.0
    if length_anomaly and repeat and array_frac >= 0.8 and ref_len and L / ref_len >= 1.6:
        out["anomaly_type"] = "concatemer"
        out["suggested_trim_region"] = f"{ref_len + 1}-{L}"
        out["curation_suggestion"] = (
            f"Assembly is {round(L / ref_len, 2)}x the reference (tandem genome "
            f"duplication / concatemer). Collapse to a single monomer (~{ref_len} bp), "
            f"removing positions {out['suggested_trim_region']} (~{excess} bp)."
        )
    elif length_anomaly and repeat and in_cr in ("yes", "NA"):
        out["anomaly_type"] = "control_region_repeat"
        # Trim the distal (lowest-coverage) end of the array by the excess length
        # to bring the assembly back toward the reference length.
        trim_bp = min(excess, repeat["array_len"] - repeat["unit_bp"])
        if trim_bp > 0:
            trim_start = repeat["array_end"] - trim_bp + 1
            out["suggested_trim_region"] = f"{trim_start}-{repeat['array_end']}"
            cr_txt = "control-region " if in_cr == "yes" else ""
            out["curation_suggestion"] = (
                f"Assembly {excess} bp longer than reference; {cr_txt}tandem repeat "
                f"(~{repeat['unit_bp']} bp unit x ~{repeat['copies']} copies, "
                f"{repeat['identity']}% id) at {out['repeat_region']}. Review for length "
                f"heteroplasmy; to approach reference length collapse distal array, "
                f"removing ~{trim_bp} bp at {out['suggested_trim_region']}."
            )
    elif length_anomaly:
        out["anomaly_type"] = "unresolved"
        out["curation_suggestion"] = (
            f"Assembly {excess} bp longer than reference with no clear tandem repeat; "
            f"review for duplication / NUMT contamination."
        )
    return out


def seq_fasta(seq, workdir):
    """Write a sequence to a FASTA in workdir and return the path (memoised)."""
    path = os.path.join(workdir, "assembly.fasta")
    if not os.path.exists(path):
        with open(path, "w") as fh:
            fh.write(">assembly\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")
    return path


def count_junction_spanning_reads(fasta_seq, bam, sample, threads, min_overhang, workdir):
    """Remap the mapped HiFi reads to a doubled reference and count reads that
    bridge the linearisation point. Returns (spanning_reads, mapped_reads)."""
    L = len(fasta_seq)
    doubled = os.path.join(workdir, "doubled.fasta")
    with open(doubled, "w") as fh:
        fh.write(f">{sample}_doubled\n")
        seq2 = fasta_seq + fasta_seq
        for i in range(0, len(seq2), 80):
            fh.write(seq2[i : i + 80] + "\n")

    reads = os.path.join(workdir, "mapped_reads.fasta")
    # Primary, mapped reads only (exclude unmapped/secondary/supplementary).
    with open(reads, "w") as out:
        run(["samtools", "fasta", "-F", "0x904", bam], stdout=out)

    if os.path.getsize(reads) == 0:
        log("No reads extracted from BAM; cannot run read-spanning test.")
        return 0, 0

    remap = os.path.join(workdir, "remap.bam")
    # minimap2 HiFi preset -> primary alignments only -> sorted BAM.
    pipeline = (
        f"minimap2 -ax map-hifi -t {threads} {doubled} {reads} "
        f"| samtools view -b -F 0x904 - "
        f"| samtools sort -@ {threads} -o {remap} -"
    )
    run(pipeline, shell=True)
    run(["samtools", "index", remap])

    spanning = set()
    mapped = set()
    view = subprocess.run(
        ["samtools", "view", "-F", "0x904", remap],
        check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,  # 'text=' alias; portable to the container's older Python
    )
    for line in view.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        qname, pos, cigar = fields[0], fields[3], fields[5]
        if cigar == "*":
            continue
        mapped.add(qname)
        ref_start = int(pos) - 1  # 0-based
        ref_end = ref_start + ref_aligned_length(cigar)  # exclusive
        # Junction at coordinate L. Require the alignment to start in the first
        # copy and to extend min_overhang bp past the junction on both sides.
        if ref_start < L and ref_start <= L - min_overhang and ref_end >= L + min_overhang:
            spanning.add(qname)
    return len(spanning), len(mapped)


def parse_contig_stats(path):
    """Parse a MitoHiFi contig-stats table.

    Returns a dict with: comments (list), header (list), rows (list of list),
    col index map, the final_mitogenome row index, the current was_circular
    value, and the list of hifiasm source-contig ids for the final mitogenome.
    """
    comments, header, rows = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                comments.append(line)
            elif header is None:
                header = line.split("\t")
            else:
                if line.strip():
                    rows.append(line.split("\t"))
    if header is None:
        raise ValueError(f"No header found in {path}")

    idx = {name: i for i, name in enumerate(header)}
    for required in ("contig_id", "was_circular"):
        if required not in idx:
            raise ValueError(f"Column '{required}' missing from {path}")

    ci, wi = idx["contig_id"], idx["was_circular"]
    ai = idx.get("annotation_file")
    final_idx = None
    final_circ = None
    source_contigs = []
    for i, row in enumerate(rows):
        if len(row) <= max(ci, wi):
            continue
        if row[ci] == "final_mitogenome":
            final_idx = i
            final_circ = row[wi]
        # Source hifiasm contig(s) annotated as the final mitogenome.
        elif ai is not None and len(row) > ai:
            ann = os.path.basename(row[ai])
            if ann == "final_mitogenome.gb":
                source_contigs.append(row[ci])
    return {
        "comments": comments,
        "header": header,
        "rows": rows,
        "idx": idx,
        "final_idx": final_idx,
        "final_circular": final_circ,
        "source_contigs": source_contigs,
    }


def coerce_bool(value):
    if value is None:
        return None
    s = str(value).strip().lower()
    if s in {"true", "t", "1", "yes", "y", "circular"}:
        return True
    if s in {"false", "f", "0", "no", "n", "linear"}:
        return False
    return None


def write_corrected_stats(stats, set_circular_for, out_path):
    """Write the stats table, forcing was_circular=True for the given contig ids."""
    ci = stats["idx"]["contig_id"]
    wi = stats["idx"]["was_circular"]
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as out:
        for c in stats["comments"]:
            out.write(c + "\n")
        out.write("\t".join(stats["header"]) + "\n")
        for row in stats["rows"]:
            if len(row) > max(ci, wi) and row[ci] in set_circular_for:
                row = list(row)
                row[wi] = "True"
            out.write("\t".join(row) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="Final mitogenome FASTA.")
    p.add_argument("--contigs-stats", required=True, help="MitoHiFi contig-stats (.with_coverage) TSV.")
    p.add_argument("--bam", help="HiFi-vs-final_mitogenome sorted BAM for the read-spanning test.")
    p.add_argument("--gb", help="Final mitogenome GenBank (for control-region location of repeats).")
    p.add_argument("--sample", required=True, help="Sample / assembly prefix.")
    p.add_argument("--out-stats", required=True, help="Output corrected contig-stats TSV.")
    p.add_argument("--out-evidence", required=True, help="Output circularity-check evidence TSV.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-spanning-reads", type=int, default=2, help="Reads bridging the junction to call circular.")
    p.add_argument("--min-overhang", type=int, default=200, help="Min bp an alignment must extend past the junction each side.")
    p.add_argument("--reference-length", type=int, help="Override related-reference length (else parsed from contig-stats header).")
    p.add_argument("--max-length-ratio", type=float, default=1.10, help="assembly/reference length ratio above which to flag a length anomaly.")
    p.add_argument("--min-repeat-length", type=int, default=200, help="Min self-alignment length (bp) to consider a tandem repeat.")
    p.add_argument("--min-repeat-identity", type=float, default=70.0, help="Min %% identity for a tandem-repeat self-alignment.")
    args = p.parse_args()

    error_note = ""
    stats = parse_contig_stats(args.contigs_stats)
    mitohifi_circ = coerce_bool(stats["final_circular"])

    # Signal 2: hifiasm c/l flag from the source contig name(s).
    source_contigs = stats["source_contigs"]
    hifiasm_circ = any(c.endswith("c") for c in source_contigs) if source_contigs else None
    rep_contig = ";".join(source_contigs) if source_contigs else "NA"

    # Signal 1: read-spanning junction test (skipped if already circular, or no BAM).
    spanning_reads = "NA"
    mapped_reads = "NA"
    read_span_circ = None
    _, seq = read_fasta(args.fasta)

    if mitohifi_circ is True:
        # MitoHiFi already called it circular; nothing to correct.
        log("MitoHiFi already reports circular; passing through unchanged.")
    elif not seq:
        error_note = "empty_fasta"
        log("Final FASTA is empty; cannot run read-spanning test.")
    elif args.bam and os.path.exists(args.bam) and os.path.getsize(args.bam) > 0:
        try:
            with tempfile.TemporaryDirectory() as workdir:
                n_span, n_map = count_junction_spanning_reads(
                    seq, args.bam, args.sample, args.threads, args.min_overhang, workdir
                )
            spanning_reads, mapped_reads = n_span, n_map
            read_span_circ = n_span >= args.min_spanning_reads
        except subprocess.CalledProcessError as e:
            error_note = "remap_failed"
            log(f"Read-spanning test failed ({e}); falling back to hifiasm flag.")
        except Exception as e:  # pragma: no cover - defensive
            error_note = "remap_error"
            log(f"Read-spanning test error ({e}); falling back to hifiasm flag.")
    else:
        error_note = "bam_missing"
        log("No usable BAM; relying on hifiasm flag only.")

    # Verdict: circular if MitoHiFi already closed it, or either independent
    # signal is positive. MitoHiFi's own verdict must be carried through here:
    # when it reports circular the read-spanning test is skipped (read_span_circ
    # stays None) and a hifiasm '...l' unitig gives hifiasm_circ=False, so
    # without this the verdict would collapse to False for a genuinely circular
    # genome and propagate a linear topology downstream (meta.circular).
    verdict = (mitohifi_circ is True) or bool(read_span_circ) or bool(hifiasm_circ)
    corrected = bool(mitohifi_circ is not True and verdict)

    if corrected:
        targets = set(source_contigs) | {"final_mitogenome"}
        write_corrected_stats(stats, targets, args.out_stats)
        log(f"CORRECTED was_circular -> True (spanning_reads={spanning_reads}, hifiasm_circular={hifiasm_circ}).")
    else:
        # Pass the table through unchanged (force no contig ids).
        write_corrected_stats(stats, set(), args.out_stats)
        log("No correction applied; stats passed through unchanged.")

    # Length / tandem-repeat assessment (runs for every assembly, regardless of
    # the circularity verdict).
    ref_len = args.reference_length or parse_reference_length(args.contigs_stats)
    lr = {}
    if seq:
        with tempfile.TemporaryDirectory() as workdir:
            lr = assess_length_and_repeats(seq, ref_len, args.gb, workdir, args)
    if lr.get("curation_suggestion", "none") != "none":
        log(f"LENGTH/REPEAT: {lr['curation_suggestion']}")

    def fmt(v):
        return "NA" if v is None else ("True" if v is True else ("False" if v is False else str(v)))

    lr_cols = [
        "assembly_length", "reference_length", "length_ratio", "excess_bp", "length_anomaly",
        "tandem_repeat", "repeat_region", "repeat_unit_bp", "repeat_copies", "repeat_identity",
        "repeat_in_control_region", "anomaly_type", "suggested_trim_region", "curation_suggestion",
    ]

    os.makedirs(os.path.dirname(args.out_evidence) or ".", exist_ok=True)
    with open(args.out_evidence, "w") as out:
        out.write(
            "sample\tsource_contig\tmitohifi_was_circular\thifiasm_contig_circular\t"
            "junction_spanning_reads\tmapped_reads\tmin_spanning_reads\tmin_overhang\t"
            "read_spanning_circular\tfinal_verdict_circular\tflag_corrected\tnote\t"
            + "\t".join(lr_cols) + "\n"
        )
        out.write(
            f"{args.sample}\t{rep_contig}\t{fmt(mitohifi_circ)}\t{fmt(hifiasm_circ)}\t"
            f"{spanning_reads}\t{mapped_reads}\t{args.min_spanning_reads}\t{args.min_overhang}\t"
            f"{fmt(read_span_circ)}\t{fmt(verdict)}\t{'yes' if corrected else 'no'}\t{error_note or 'ok'}\t"
            + "\t".join(str(lr.get(c, "NA")) for c in lr_cols) + "\n"
        )

    log("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
