#!/usr/bin/env python3
"""Diff which reference SELECT_REFERENCE_DB picks, between two builds of a group DB.

WHY THIS EXISTS
---------------
Rebuilding a group database (lifting `refseq[filter]`, refreshing from NCBI, changing
the per-organism cap) changes which record wins the bitscore ranking in
select_reference_db.py, and therefore the annotation reference, the QC baseline and the
top-n GetOrganelle reseed seed. The record count alone does not say whether that is an
improvement: a database can triple in size and move no pick, or add one congeneric
record and rescue a sample that had nothing closer than another order.

This measures it directly. For every assembly in a corpus it runs the real
select_reference_db.py twice -- once against the old database directory, once against
the new one -- and grades both winners against the sample's own taxonomy. It reuses the
pipeline's own code rather than reimplementing the ranking, so what it reports is what
the pipeline would do.

It is an audit tool, not a pipeline step, and it is deliberately read-only: it writes
one TSV and a scratch working directory, and never touches assets/refdb.

WHAT A VERDICT MEANS
--------------------
    equal_record  same accession won both times -- the rebuild did not move this sample
    same          different record, same taxonomic tier (usually a sibling species)
    closer        the new database has a taxonomically nearer reference. The point.
    further       the new winner is taxonomically MORE distant than the old one

`further` is worth reading before committing a rebuild. It usually means a longer or
better-covering record outscored a nearer one on total bitscore, which is a statement
about the ranking rather than about the database -- but it can also mean a new record
is mis-annotated in GenBank, so look at the row rather than at the count.

Grading is by taxonomy, not sequence, and is deliberately invertebrate-shaped:
reference_divergence_check.order_from_lineage() only recognises vertebrate orders (the
'iformes' suffix), which never matches Anomura, Platyctenida or Scleractinia. Here the
sample's own class/order/family from the samplesheet are matched against the reference
record's NCBI lineage instead, so the order tier works for invertebrates.

WHERE IT RUNS
-------------
select_reference_db.py needs blastn and biopython, which are not on a Setonix login
node. Run this inside the MITOS container, the same one the module uses:

    singularity exec --env PYTHONPATH=/path/to/repo/bin \\
        $SING/depot.galaxyproject.org-singularity-mitos-2.1.10--pyhdfd78af_0.img \\
        python3 bin/audit_reference_selection_diff.py ...

Usage:
    audit_reference_selection_diff.py \\
        --samplesheet /scratch/.../invert_test_panel/panel_all.csv \\
        --corpus /scratch/.../invert_test_panel/out_panel_all/mitogenomes \\
        --old-refdb-root assets/refdb \\
        --new-refdb-root /scratch/.../refdb_widened \\
        --workdir /scratch/.../refdb_widened/_diff_work \\
        --out /scratch/.../refdb_widened/selection_diff.tsv
"""
import argparse
import collections
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parent
REPO_ROOT = BIN_DIR.parent
sys.path.insert(0, str(BIN_DIR))

# Reused rather than re-typed: parse_groovy_class_groups() reads the class -> group
# map straight out of InvertTaxonGroups.groovy, which is the same stage-1 narrowing
# the pipeline does, and parse_reference() is the reference GenBank reader the
# divergence check already uses.
from build_origin_anchor_table import parse_groovy_class_groups  # noqa: E402
from reference_divergence_check import parse_reference  # noqa: E402
from species_name_utils import genus_of  # noqa: E402

GROOVY_SETS = REPO_ROOT / "lib" / "InvertTaxonGroups.groovy"

# Taxonomic tiers, nearest first. The integer is what makes closer/further orderable;
# UNKNOWN is negative so "we could not grade it" never reads as an improvement over
# a graded result. DISTANT is the floor rather than "from another phylum": a group
# database only ever holds one phylum, so a DISTANT pick is a record that shares the
# phylum and nothing below it -- Tjalfiella (Tentaculata) seeded from Beroe (Nuda) is
# the case this whole audit was written to catch.
TIERS = {"CONGENERIC": 4, "CONFAMILIAL": 3, "SAME_ORDER": 2,
         "SAME_CLASS": 1, "DISTANT": 0, "UNKNOWN": -1}

COLUMNS = ["sample", "group", "assembly",
           "old_state", "old_acc", "old_organism", "old_family", "old_tier",
           "old_cov", "old_pid",
           "new_state", "new_acc", "new_organism", "new_family", "new_tier",
           "new_cov", "new_pid",
           "verdict", "note"]


def sample_genus(nominal, family):
    """The sample's genus, or '' when the label is not a genus-level name.

    OceanOmics invert labels are often coarse -- INV01's nominal_species_id is
    'Acanthogorgiidae', a family. Treating that as a genus would compare a family name
    against a reference genus and grade every reference DISTANT, so a label that is
    the family (or just looks like one) yields no genus and the grading falls through
    to the family tier, which is the finest rank the label actually asserts.
    """
    g = genus_of(nominal or "")
    if not g:
        return ""
    if g.lower().endswith("idae") or (family and g.lower() == family.lower()):
        return ""
    return g


def grade(sample_row, ref_organism, ref_lineage):
    """(tier, score) for a reference against a sample, from taxonomy alone.

    Pure, so it is unit-testable without biopython or a database. Matches the sample's
    class/order/family against the reference record's NCBI lineage, which is the part
    reference_divergence_check cannot do for invertebrates.
    """
    ref_genus = genus_of(ref_organism or "")
    lineage = {t.strip().lower() for t in (ref_lineage or []) if t.strip()}
    fam = (sample_row.get("family") or "").strip()
    order = (sample_row.get("order") or "").strip()
    cls = (sample_row.get("class") or "").strip()
    genus = sample_genus(sample_row.get("nominal_species_id"), fam)

    if genus and ref_genus and genus.lower() == ref_genus.lower():
        tier = "CONGENERIC"
    elif fam and fam.lower() in lineage:
        tier = "CONFAMILIAL"
    elif order and order.lower() in lineage:
        tier = "SAME_ORDER"
    elif cls and cls.lower() in lineage:
        tier = "SAME_CLASS"
    elif ref_genus:
        tier = "DISTANT"
    else:
        tier = "UNKNOWN"
    return tier, TIERS[tier]


def selected(result):
    """Did select_reference_db.py actually choose a record? Covers SELECTED and
    SELECTED_LOW_CONFIDENCE, and excludes NONE, MISSING and the error states."""
    return result.get("state", "").startswith("SELECTED")


def verdict_of(old, new):
    """closer / further / same / equal_record / no_selection, from two results.

    'no_selection' is separate from 'same' on purpose. A sample that aligns to nothing
    in either database (INV08_TJALFIELLA: 'no DB record aligned to the assembly', a
    4.8 kb fragment at ~4x) is not evidence that the rebuild changed nothing -- it is
    the absence of a measurement, and counting it as 'same' would quietly pad the
    unchanged column with samples the audit could not speak to at all.
    """
    if not selected(old) and not selected(new):
        return "no_selection"
    if old["acc"] and old["acc"] == new["acc"]:
        return "equal_record"
    if new["score"] > old["score"]:
        return "closer"
    if new["score"] < old["score"]:
        return "further"
    return "same"


def read_samplesheets(paths):
    """sample id -> taxonomy row, from the pipeline's own samplesheet CSV."""
    rows = {}
    for p in paths:
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh):
                sid = (row.get("sample") or "").strip()
                if sid:
                    rows[sid] = row
    if not rows:
        sys.exit(f"[refdiff] no samples parsed from {', '.join(map(str, paths))}")
    return rows


def find_assemblies(corpus_dirs):
    """sample id -> the largest published mitogenome FASTA for that sample.

    A published sample directory holds one subdirectory per assembly prefix
    (<prefix>/mtdna/<prefix>.fasta), and a reseeded sample has both its empty first
    pass and its reseed. The largest non-empty file is the assembly the run actually
    produced, and it is what SELECT_REFERENCE_DB would be ranking against. Empty
    files are reported as skipped rather than silently dropped -- an empty first pass
    is exactly the case select_fallback_seed exists for, and it has no assembly to
    BLAST, so no reference diff is possible for it.
    """
    best, empty = {}, {}
    for root in corpus_dirs:
        root = Path(root)
        for fa in sorted(root.glob("*/*/mtdna/*.fasta")):
            sid = fa.relative_to(root).parts[0]
            size = fa.stat().st_size
            if size == 0:
                empty.setdefault(sid, fa)
                continue
            if sid not in best or size > best[sid].stat().st_size:
                best[sid] = fa
    return best, {s: f for s, f in empty.items() if s not in best}


def parse_status(path):
    """(state, detail) from a select_reference_db.py status line."""
    if not path.exists():
        return "MISSING", ""
    text = path.read_text().strip()
    if not text:
        return "MISSING", ""
    state, _, detail = text.partition("\t")
    return state.strip(), detail.strip()


def parse_metrics(detail):
    """cov and pid out of a 'SELECTED' detail line, as strings ('' when absent)."""
    out = {"cov": "", "pid": ""}
    for token in detail.split():
        for key in ("cov", "pid"):
            if token.startswith(f"{key}="):
                out[key] = token.split("=", 1)[1]
    return out


def select(assembly, refdb_dir, group, workdir, tag, extra_args):
    """Run select_reference_db.py in reference mode and read back what it chose."""
    out_gb = workdir / f"{tag}.reference.gb"
    out_status = workdir / f"{tag}.reference_select.txt"
    cmd = [sys.executable, str(BIN_DIR / "select_reference_db.py"),
           "--assembly", str(assembly),
           "--refdb-dir", str(refdb_dir),
           "--group", group,
           "--out-gb", str(out_gb),
           "--out-status", str(out_status), *extra_args]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    state, detail = parse_status(out_status)
    result = {"state": state, "acc": "", "organism": "", "family": "",
              "tier": "", "score": TIERS["UNKNOWN"], "lineage": [],
              **parse_metrics(detail)}
    if proc.returncode != 0:
        # select_reference_db.py always exits 0 by design, so a non-zero code is an
        # environment problem (no blastn, unreadable database) and must not be
        # graded as if it were a taxonomic result.
        result["state"] = f"ERROR({proc.returncode})"
        result["note"] = (proc.stderr or "").strip().splitlines()[-1:] or [""]
        return result
    if out_gb.exists() and out_gb.stat().st_size > 0:
        try:
            organism, _genus, family, _order, lineage = parse_reference(out_gb)
            result.update(organism=organism, family=family, lineage=lineage)
        except Exception as exc:
            result["state"] = f"UNREADABLE_GB({exc})"
    # 'SELECTED <acc> <organism> [<family>] cov=.. pid=..'
    result["acc"] = detail.split()[0] if detail and state.startswith("SELECTED") else ""
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", type=Path, action="append", required=True,
                    help="Pipeline samplesheet CSV giving each sample's "
                         "nominal_species_id/class/order/family. Repeatable.")
    ap.add_argument("--corpus", type=Path, action="append", required=True,
                    help="A published mitogenomes/ directory to take assemblies from. "
                         "Repeatable.")
    ap.add_argument("--old-refdb-root", type=Path,
                    default=REPO_ROOT / "assets" / "refdb",
                    help="Root holding the CURRENT <group>/ database directories.")
    ap.add_argument("--new-refdb-root", type=Path, required=True,
                    help="Root holding the REBUILT <group>/ database directories.")
    ap.add_argument("--group", action="append",
                    help="Only audit these groups. Repeatable. Default: all.")
    # Not given defaults here on purpose. These only decide whether the winner is
    # labelled SELECTED or SELECTED_LOW_CONFIDENCE, never which record wins, and
    # conf/modules.config leaves them unset for the annotation invocation -- so the
    # audit inherits select_reference_db.py's own values rather than quietly
    # grading against a different confidence bar than the pipeline uses.
    ap.add_argument("--min-cov", type=float,
                    help="Override select_reference_db.py's confidence threshold.")
    ap.add_argument("--min-pid", type=float,
                    help="Override select_reference_db.py's confidence threshold.")
    ap.add_argument("--workdir", type=Path,
                    help="Where the per-sample GenBank/status files are written. "
                         "Default: a temporary directory, removed on exit.")
    ap.add_argument("--out", type=Path, required=True, help="Output TSV.")
    args = ap.parse_args()

    extra_args = []
    if args.min_cov is not None:
        extra_args += ["--min-cov", str(args.min_cov)]
    if args.min_pid is not None:
        extra_args += ["--min-pid", str(args.min_pid)]

    class_groups = parse_groovy_class_groups(GROOVY_SETS)
    samples = read_samplesheets(args.samplesheet)
    assemblies, empties = find_assemblies(args.corpus)
    wanted = set(args.group) if args.group else None

    tmp = None
    if args.workdir:
        workdir = args.workdir
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        tmp = tempfile.mkdtemp(prefix="refdiff_")
        workdir = Path(tmp)

    rows, tally, skipped = [], collections.Counter(), collections.Counter()
    try:
        for sid in sorted(set(samples) | set(assemblies)):
            meta = samples.get(sid)
            if meta is None:
                skipped["no_samplesheet_row"] += 1
                continue
            group = class_groups.get((meta.get("class") or "").strip().lower())
            if group is None:
                skipped["class_in_no_group"] += 1
                continue
            if wanted and group not in wanted:
                continue
            assembly = assemblies.get(sid)
            if assembly is None:
                skipped["empty_assembly" if sid in empties else "no_assembly"] += 1
                continue
            old_dir = args.old_refdb_root / group
            new_dir = args.new_refdb_root / group
            if not old_dir.is_dir() or not new_dir.is_dir():
                skipped["database_missing"] += 1
                continue

            old = select(assembly, old_dir, group, workdir, f"{sid}.old", extra_args)
            new = select(assembly, new_dir, group, workdir, f"{sid}.new", extra_args)
            for res in (old, new):
                res["tier"], res["score"] = grade(meta, res["organism"], res["lineage"])
            verdict = verdict_of(old, new)
            tally[verdict] += 1
            rows.append({
                "sample": sid, "group": group, "assembly": str(assembly),
                "old_state": old["state"], "old_acc": old["acc"],
                "old_organism": old["organism"], "old_family": old["family"],
                "old_tier": old["tier"], "old_cov": old["cov"], "old_pid": old["pid"],
                "new_state": new["state"], "new_acc": new["acc"],
                "new_organism": new["organism"], "new_family": new["family"],
                "new_tier": new["tier"], "new_cov": new["cov"], "new_pid": new["pid"],
                "verdict": verdict,
                "note": "" if selected(old) and selected(new) else "not_selected",
            })
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"[refdiff] {len(rows)} samples audited -> {args.out}")
    for verdict in ("closer", "same", "equal_record", "further", "no_selection"):
        if tally[verdict]:
            print(f"[refdiff]   {verdict:12s} {tally[verdict]}")
    for reason, n in sorted(skipped.items()):
        print(f"[refdiff]   skipped {reason}: {n}", file=sys.stderr)
    # Every 'closer' row is a result worth quoting in the CHANGELOG, and every
    # 'further' row is one to read before committing the rebuild, so name them here
    # rather than making the reader grep the TSV.
    for row in rows:
        if row["verdict"] in ("closer", "further"):
            print(f"[refdiff]   {row['verdict']:8s} {row['sample']:22s} "
                  f"{row['old_acc']} {row['old_organism']} [{row['old_tier']}] -> "
                  f"{row['new_acc']} {row['new_organism']} [{row['new_tier']}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
