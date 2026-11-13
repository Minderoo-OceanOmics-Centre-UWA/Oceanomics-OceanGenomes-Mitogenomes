#!/usr/bin/env python3
"""Pick a sample's reference records from a curated group DB by SEQUENCE.

This is the second of the two narrowing stages that resolve a seed/reference.
Stage 1 picks the group (InvertTaxonGroups.seedDbGroup: class -> one of the
databases under assets/refdb/); this stage picks records out of that group by
BLASTing the sample's own assembly against it. It is label-free, so a wrong or
coarse species label can no longer decide the answer.

Two output modes, both driven by the same ranking:

  seed mode (--out-seed-fasta / --out-label-fasta, optionally --out-gb)
      The top --top-n records become GETORGANELLE_RESEED's seed (-s) and their
      genes its label database (--genes). Without this the reseed was handed the
      WHOLE group -- 221 anthozoan genomes, 850 molluscan -- which recruits reads
      from across the phylum, halves effective coverage and shatters the graph:
      the case that motivated this was a coral whose reseed went from 2 scaffolds
      to 12. Top-n rather than the single best because from a small fragmented
      first pass the pick is reliable at order/subclass level, not species.

  reference mode (--out-gb)
      The single best record becomes the annotation reference for
      CORAL_ANNOTATION_FIX, and the reference GETORGANELLE_CHECK and
      REFERENCE_RELEVANCE grade the assembly against.

Ranking = total bitscore against the assembly (closest + most-covering). The
winner is reported with its assembly coverage and coverage-weighted identity;
below --min-cov/--min-pid it is still emitted but marked LOW_CONFIDENCE so the
downstream relevance check can surface it (e.g. a lineage with no DB record).

The DB is addressed as --refdb-dir + --group, reading the tracked
<group>_mito_refdb.{fasta,label.fasta,manifest.tsv,features.tsv}; the reference
record is rebuilt from those by refdb_record.py rather than sliced out of a
stored .gb, which is what lets every group have a reference and not just the one
group whose .gb was small enough to track. --db-gb remains supported for a
self-contained multi-record GenBank (the module tests use one).

Always exits 0. Runs in the MITOS2 BioContainer (blastn + biopython).

Usage:
    select_reference_db.py --assembly asm.fa --refdb-dir assets/refdb/anthozoa \
        --group anthozoa --out-status PREFIX.reference_select.txt \
        --out-seed-fasta PREFIX.seed.fasta --out-label-fasta PREFIX.genedb.fasta
    select_reference_db.py --assembly asm.fa --db-gb mini_db.gb \
        --out-gb PREFIX.reference.gb --out-status PREFIX.reference_select.txt
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

# bin/ is on PATH at runtime, so a sibling import resolves via sys.path[0]
# (same pattern as bin/process_files.py -> orf_utils).
import refdb_record
from refdb_record import subset_label_db, union_coverage


def family_of(rec):
    return next((t for t in rec.annotations.get("taxonomy", []) or [] if t.endswith("idae")), "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True, type=Path)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--refdb-dir", type=Path,
                     help="assets/refdb/<group>/ holding the tracked artifacts")
    src.add_argument("--db-gb", type=Path, help="self-contained multi-record GenBank")
    ap.add_argument("--group", help="group name (required with --refdb-dir)")
    ap.add_argument("--out-status", required=True, type=Path)
    ap.add_argument("--out-gb", type=Path, help="reference mode: single best record")
    ap.add_argument("--out-seed-fasta", type=Path, help="seed mode: top-n genomes")
    ap.add_argument("--out-label-fasta", type=Path, help="seed mode: their label records")
    ap.add_argument("--top-n", type=int, default=5)
    ap.add_argument("--min-cov", type=float, default=0.60)
    ap.add_argument("--min-pid", type=float, default=88.0)
    args = ap.parse_args()

    if args.refdb_dir and not args.group:
        ap.error("--group is required with --refdb-dir")
    seed_mode = bool(args.out_seed_fasta or args.out_label_fasta)
    if seed_mode and not (args.out_seed_fasta and args.out_label_fasta):
        ap.error("seed mode needs both --out-seed-fasta and --out-label-fasta")

    def status(state, msg):
        args.out_status.write_text(f"{state}\t{msg}\n")
        print(f"[select_reference_db] {state}: {msg}", file=sys.stderr)
        sys.exit(0)

    asm_recs = list(SeqIO.parse(str(args.assembly), "fasta"))
    if not asm_recs:
        status("NONE", f"empty assembly {args.assembly.name}")
    asm_len = sum(len(r.seq) for r in asm_recs)

    # BLAST subject: the group's tracked genome FASTA, or the records of a
    # self-contained GenBank written out to one.
    tmp = tempfile.TemporaryDirectory()
    if args.refdb_dir:
        db_fa, manifest_p, features_p = refdb_record.refdb_paths(args.refdb_dir, args.group)
        if not db_fa.exists():
            status("NONE", f"no database at {db_fa}")
        db_label = Path(str(db_fa).replace(".fasta", ".label.fasta"))
    else:
        gb_recs = list(SeqIO.parse(str(args.db_gb), "genbank"))
        if not gb_recs:
            status("NONE", f"empty DB {args.db_gb.name}")
        db_fa = Path(tmp.name) / "db.fasta"
        SeqIO.write(gb_recs, str(db_fa), "fasta")
        manifest_p = features_p = db_label = None

    out = subprocess.run(
        ["blastn", "-query", str(args.assembly), "-subject", str(db_fa),
         "-evalue", "1e-10", "-max_hsps", "20",
         "-outfmt", "6 sseqid qstart qend pident length bitscore"],
        capture_output=True, text=True).stdout

    # Aggregate per DB subject: total bitscore (ranking), assembly coverage +
    # weighted identity (the guard on the winner).
    agg = {}
    for line in out.strip().splitlines():
        sid, qs, qe, pid, ln, bits = line.split("\t")
        qs, qe, pid, ln, bits = int(qs), int(qe), float(pid), int(ln), float(bits)
        a = agg.setdefault(sid, {"bits": 0.0, "iv": [], "alnbp": 0, "wpid": 0.0})
        a["bits"] += bits
        a["iv"].append((min(qs, qe), max(qs, qe)))
        a["alnbp"] += ln
        a["wpid"] += ln * pid
    if not agg:
        status("NONE", "no DB record aligned to the assembly")

    ranked = sorted(agg, key=lambda k: agg[k]["bits"], reverse=True)
    best_id = ranked[0]
    a = agg[best_id]
    cov = union_coverage(a["iv"]) / asm_len
    pid = a["wpid"] / a["alnbp"] if a["alnbp"] else 0.0

    # --- seed mode: top-n genomes + the matching slice of the label database ---
    if seed_mode:
        if db_label is None or not db_label.exists():
            status("NONE", "seed mode needs --refdb-dir (no label database available)")
        chosen = ranked[:max(1, args.top_n)]
        picked = [r for r in SeqIO.parse(str(db_fa), "fasta") if r.id in set(chosen)]
        SeqIO.write(picked, str(args.out_seed_fasta), "fasta")
        n_label = subset_label_db(db_label, chosen, args.out_label_fasta)
        if not picked or not n_label:
            status("NONE", f"no seed records materialised for {','.join(chosen)}")
        fams = {}
        if manifest_p and manifest_p.exists():
            rows = refdb_record.read_manifest(manifest_p)
            fams = {c: rows.get(c, {}).get("family", "") for c in chosen}
        # Emit the single best record too when asked. The reseed wants both: the
        # top-n as its seed, and the best one as the reference GETORGANELLE_CHECK
        # and REFERENCE_RELEVANCE grade the assembly against -- which invertebrates
        # never had, so every invert check recorded note=no_reference with NA
        # coverage. One ranking already answers both, so this costs no second BLAST.
        # Unreachable in practice -- best_id came from BLASTing against this very
        # FASTA, so the record is there -- but a SELECTED_SEED status with no
        # reference beside it would leave the sample in neither the reference nor
        # the placeholder channel, and so drop it from the run. Report NONE instead
        # and let it fall back to its first pass.
        if args.out_gb:
            rec = refdb_record.materialize(args.refdb_dir, args.group, best_id)
            if rec is None:
                status("NONE", f"{best_id} ranked first but could not be materialised "
                               f"from {args.refdb_dir}")
            SeqIO.write([rec], str(args.out_gb), "genbank")

        detail = ", ".join(f"{c}[{fams.get(c, '')}]" for c in chosen)
        status("SELECTED_SEED",
               f"{len(picked)} of {len(agg)} aligned records, {n_label} label seqs: "
               f"{detail} best_cov={cov:.2f} best_pid={pid:.1f}")

    # --- reference mode: the single best record -------------------------------
    if args.refdb_dir:
        rec = refdb_record.materialize(args.refdb_dir, args.group, best_id)
        if rec is None:
            status("NONE", f"{best_id} not found in {args.refdb_dir}")
    else:
        rec = next(r for r in gb_recs if r.id == best_id)

    SeqIO.write([rec], str(args.out_gb), "genbank")
    organism = rec.annotations.get("organism", best_id)
    fam = family_of(rec)
    confident = cov >= args.min_cov and pid >= args.min_pid
    state = "SELECTED" if confident else "SELECTED_LOW_CONFIDENCE"
    status(state, f"{best_id} {organism} [{fam}] cov={cov:.2f} pid={pid:.1f}")


if __name__ == "__main__":
    main()
