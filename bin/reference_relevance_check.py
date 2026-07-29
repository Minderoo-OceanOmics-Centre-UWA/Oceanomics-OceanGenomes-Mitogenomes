#!/usr/bin/env python3
"""Grade how well the mitogenome reference chosen for a sample corresponds to its
assembly.

The reference is resolved by MITOHIFI_FINDMITOREFERENCE from the sample's species
*label* (reference_species_id / nominal_species_id). A wrong or coarse label
therefore yields a wrong-family reference, which silently degrades both
GetOrganelle seeding and the anthozoan annotation fix (e.g. a Merulinid sample
labelled "Montipora grisea" gets an Acroporidae reference its 16S/nad5 cannot be
transferred from). MITOS/coral_fix then just emit a partial annotation with no
obvious cause.

This check is label-free and taxonomy-DB-free: it BLASTs the chosen reference
genome against the assembled mitogenome and measures how much of the reference is
present in the assembly, and at what identity.

Three states, because "not a close relative" and "not the right molecule" are
different problems with different remedies:

    PASS       reference covers the assembly at the expected identity
    DIVERGENT  reference corresponds to the assembly but is distant -- the
               assembly is fine; a closer reference would seed/annotate better
    MISMATCH   reference neither covers nor matches: almost certainly the wrong
               reference for this sample (bad species label)

Three things this deliberately does NOT do, each learned from real data:

* **Coverage is normalised by REFERENCE length, not assembly length.** Assembly-
  normalised coverage conflates "bad reference" with "inflated or fragmented
  assembly": a concatemer or a control-region tandem repeat makes the assembly
  longer than any reference can cover, and a fragmented assembly pads it with
  non-mito contigs. OG778 scored 0.30 assembly-coverage against its own species'
  reference at 100% identity purely because the assembly was fragmented.
  "How much of the reference is present" is what actually measures relevance.

* **dc-megablast, not blastn's default megablast.** The reference is a related
  species at ~70-85% identity and megablast's long exact seeds miss most of those
  HSPs. Same lesson, same fix as bin/check_getorganelle.py: on OG56 megablast
  reported 0.35 coverage where dc-megablast reports 0.76.

* **The identity floor is taxon-dependent.** The original 88.0 came from coral
  data (same genus ~99.6%, same family ~96.9%, wrong family ~81.5%), but anthozoan
  mtDNA evolves far more slowly than teleost mtDNA, where *congeneric* references
  routinely sit at 78-88%. Applying the coral number to fish flagged 27 assemblies
  in one run, 10 of them against a same-genus reference. --min-pid therefore
  defaults to the vertebrate value and the module raises it for invertebrates.

A congeneric reference is never called MISMATCH (capped at DIVERGENT): the sample
and the reference are in the same genus, so the reference is the best obtainable
and low identity is biology rather than a labelling error.

Emits one line:  PASS|DIVERGENT|MISMATCH|UNKNOWN \t <details>  and always exits 0,
so a missing/odd reference can never break the run -- it just records a flag.

Usage:
    reference_relevance_check.py --assembly asm.fa --reference-gb ref.gb \
        --out PREFIX.reference_relevance.txt [--sample-species "Genus species"] \
        [--min-cov 0.70] [--min-pid 82.0]
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

# Shared with reference_divergence_check.py so the BLAST grader and the taxonomy
# grader decide "congeneric" identically. Dependency-free, so importing it does not
# cost classify_relevance its biopython-free testability.
from species_name_utils import same_genus

# Per-HSP floors, identical to bin/check_getorganelle.py: a related-species
# reference aligns at ~70-85%, so a higher identity floor discards genuine
# coverage, and short HSPs are repeat noise rather than orthology.
MIN_HSP_LENGTH = 100
MIN_HSP_IDENTITY = 70.0


def ref_to_fasta(gb_path, fa_path):
    """Write the reference GenBank sequence(s) to FASTA; return (organism, family, length)."""
    # Imported lazily, like reference_divergence_check.parse_reference, so the pure
    # grading logic (classify_relevance) stays importable and unit-testable without
    # biopython/numpy installed.
    from Bio import SeqIO

    recs = list(SeqIO.parse(str(gb_path), "genbank"))
    if not recs:
        raise ValueError("no records in reference GenBank")
    SeqIO.write(recs, str(fa_path), "fasta")
    rec = recs[0]
    organism = rec.annotations.get("organism", rec.id)
    lineage = rec.annotations.get("taxonomy", []) or []
    family = next((t for t in lineage if t.endswith("idae")), "")
    return organism, family, sum(len(r.seq) for r in recs)


def assembly_length(fa_path):
    from Bio import SeqIO

    return sum(len(r.seq) for r in SeqIO.parse(str(fa_path), "fasta"))


def union_coverage(intervals):
    """Total length covered by a set of (lo, hi) inclusive intervals, no double count."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    covered = 0
    clo, chi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo > chi + 1:
            covered += chi - clo + 1
            clo, chi = lo, hi
        else:
            chi = max(chi, hi)
    covered += chi - clo + 1
    return covered


def classify_relevance(ref_cov, mean_pid, aln_bp, congeneric, min_cov, min_pid):
    """Grade the reference from its alignment to the assembly. Pure (no I/O), so it
    is directly unit-testable without biopython or BLAST. Returns (state, reason).

    MISMATCH requires BOTH signals to fail: a reference that covers the molecule
    but at low identity is a distant relative (DIVERGENT), and one that matches at
    high identity over only part of it is usually a partial/fragmented assembly,
    not a wrong reference. Only when it neither covers nor matches is it the wrong
    reference.
    """
    if aln_bp == 0:
        return ("DIVERGENT" if congeneric else "MISMATCH",
                "reference does not align to assembly")
    low_cov = ref_cov < min_cov
    low_pid = mean_pid < min_pid
    if low_cov and low_pid:
        if congeneric:
            return "DIVERGENT", "congeneric reference, poorly covered and low identity"
        return "MISMATCH", "reference neither covers nor matches the assembly"
    if low_cov:
        return "DIVERGENT", "reference only partly present in the assembly"
    if low_pid:
        return "DIVERGENT", "reference covers the assembly but is distantly related"
    return "PASS", "reference covers the assembly at the expected identity"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True, type=Path)
    ap.add_argument("--reference-gb", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sample-species", default="",
                    help="Sample nominal species name. When it is congeneric with the "
                         "reference organism the verdict is capped at DIVERGENT, since "
                         "a same-genus reference is the best obtainable.")
    ap.add_argument("--min-cov", type=float, default=0.70,
                    help="Flag when less than this fraction of the REFERENCE is present "
                         "in the assembly (default 0.70).")
    ap.add_argument("--min-pid", type=float, default=82.0,
                    help="Flag when the coverage-weighted mean identity is below this. "
                         "Default 82.0 suits vertebrates; the module raises it to 88.0 "
                         "for invertebrates, whose mtDNA evolves far more slowly.")
    args = ap.parse_args()

    def finish(state, msg):
        args.out.write_text(f"{state}\t{msg}\n")
        print(f"[reference_relevance] {state}: {msg}", file=sys.stderr)
        sys.exit(0)

    if not args.assembly.exists() or args.assembly.stat().st_size == 0:
        finish("UNKNOWN", f"assembly missing/empty: {args.assembly.name}")
    if not args.reference_gb.exists() or args.reference_gb.stat().st_size == 0:
        finish("UNKNOWN", f"reference missing/empty: {args.reference_gb.name}")

    try:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            ref_fa = fh.name
        organism, family, ref_len = ref_to_fasta(args.reference_gb, ref_fa)
    except Exception as exc:
        finish("UNKNOWN", f"could not parse reference {args.reference_gb.name}: {exc}")

    asm_len = assembly_length(args.assembly)
    if asm_len == 0:
        Path(ref_fa).unlink(missing_ok=True)
        finish("UNKNOWN", "assembly length 0")
    if ref_len == 0:
        Path(ref_fa).unlink(missing_ok=True)
        finish("UNKNOWN", "reference length 0")

    # qstart/qend are reference coordinates (the reference is the query), so the
    # query intervals give "how much of the reference is present in the assembly".
    out = subprocess.run(
        ["blastn", "-task", "dc-megablast", "-query", ref_fa, "-subject", str(args.assembly),
         "-evalue", "1e-5", "-outfmt", "6 sstart send length pident qstart qend"],
        capture_output=True, text=True).stdout
    Path(ref_fa).unlink(missing_ok=True)

    asm_iv, ref_iv, aln_bp, wpid = [], [], 0, 0.0
    for line in out.strip().splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        ss, se, ln, pid, qs, qe = (int(fields[0]), int(fields[1]), int(fields[2]),
                                   float(fields[3]), int(fields[4]), int(fields[5]))
        if ln < MIN_HSP_LENGTH or pid < MIN_HSP_IDENTITY:
            continue
        asm_iv.append((min(ss, se), max(ss, se)))
        ref_iv.append((min(qs, qe), max(qs, qe)))
        aln_bp += ln
        wpid += ln * pid

    ref_cov = union_coverage(ref_iv) / ref_len
    asm_cov = union_coverage(asm_iv) / asm_len
    mean_pid = (wpid / aln_bp) if aln_bp else 0.0
    congeneric = same_genus(args.sample_species, organism)

    state, reason = classify_relevance(ref_cov, mean_pid, aln_bp, congeneric,
                                       args.min_cov, args.min_pid)

    ref_desc = organism + (f" [{family}]" if family else "")
    detail = (f"{reason}; refcov={ref_cov:.2f} asmcov={asm_cov:.2f} pid={mean_pid:.1f} "
              f"ref={ref_desc}{' (congeneric)' if congeneric else ''} "
              f"thresholds=refcov>={args.min_cov:.2f},pid>={args.min_pid:.1f}")
    finish(state, detail)


if __name__ == "__main__":
    main()
