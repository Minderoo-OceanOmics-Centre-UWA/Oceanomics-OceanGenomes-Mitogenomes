#!/usr/bin/env python3
"""Re-annotate published corals against a rebuilt reference database and diff the result.

WHY THIS EXISTS
---------------
CORAL_ANNOTATION_FIX is the one consumer that copies CONTENT out of the reference rather
than merely grading against it: coral_fix_bed.py BLAST-transfers the 16S rRNA and the
intron-split nad5 (and the cox1 exons when cox1 is intron-split too) into the MITOS BED,
and mitos_to_emma.py then rebuilds the joins, the translations and the published
re-origin. So rebuilding the anthozoa database can change DEPOSITED SEQUENCE, not just a
status line, and "the new reference is in the same taxonomic tier" does not establish
that it does not.

audit_reference_selection_diff.py answers which record wins. This answers what the
annotation does with it.

WHY IT IS CHEAP
---------------
The reference enters the annotation path at exactly one point, coral_fix_bed.py --ref-gb.
MITOS2 never sees it, and MITOS2's raw output is already published under
annotation/mitos_raw/result.bed. So no MITOS2 run is needed: re-annotation is
coral_fix_bed.py followed by mitos_to_emma.py, seconds per sample, from artifacts already
on disk.

THE ONE INPUT THAT IS NOT PUBLISHED
-----------------------------------
The cox1-rotated genome the BED coordinates live in. annotation/<prefix>.fa LOOKS like it
but is the PUBLISHED RE-ORIGINED genome, a different frame (for a scleractinian it starts
at trnM, while the BED puts cox1 at 0). Annotating in the wrong frame would produce
confident nonsense, so this tool regenerates the right frame with the pipeline's own
rotate_to_cox1.py and then VERIFIES it against the BED before using it (see
check_frame()). A sample whose frame cannot be verified is reported failed, never
annotated anyway.

THREE ARMS
----------
    PUBLISHED  the annotation on disk        (as-run code, old reference)
    OLDREF     re-run                        (today's code, old reference)
    NEWREF     re-run                        (today's code, new reference)

The headline verdict is PUBLISHED vs NEWREF: what would actually change on disk. But the
published annotation was produced by older code on this branch, so that comparison alone
cannot say whether the reference or the code moved it. OLDREF costs one extra run of a
fast script and settles it:

    differs in PUBLISHED vs OLDREF  -> code drift
    differs in OLDREF   vs NEWREF  -> reference

Corals whose reference did not change are a free control: their OLDREF and NEWREF arms
must come back identical.

WHERE IT RUNS
-------------
coral_fix_bed.py and rotate_to_cox1.py need blastn/tblastn and biopython, which are not on
a Setonix login node. Run inside the MITOS container, the same one the module uses:

    singularity exec -B /scratch/pawsey1348/tpeirce -B /software/projects/pawsey1348/tpeirce \\
        $SING/depot.galaxyproject.org-singularity-mitos-2.1.10--pyhdfd78af_0.img \\
        python3 bin/audit_coral_annotation_diff.py ...

Usage:
    audit_coral_annotation_diff.py \\
        --samplesheet /scratch/.../batch-20/samplesheet/samplesheet.csv \\
        --corpus /scratch/.../batch-20/mitogenomes \\
        --refdb-dir assets/refdb/anthozoa \\
        --workdir /scratch/.../coral_annot_diff/work \\
        --out /scratch/.../coral_annot_diff/coral_annotation_diff.tsv
"""
import argparse
import collections
import csv
import hashlib
import re
import subprocess
import sys
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parent
REPO_ROOT = BIN_DIR.parent
sys.path.insert(0, str(BIN_DIR))

# The origin anchor is resolved the way the annotation subworkflow resolves it; hardcoding
# it would make the re-run diverge from the pipeline for reasons unrelated to the
# reference.
#
# The genetic code is NOT taken from create_samplesheet.resolve_genetic_code(), even
# though that is the pipeline's own resolver: create_samplesheet.py imports psycopg2 at
# module level, which is absent from both the MITOS container and refdb_venv, so it cannot
# be imported as a library here. The published annotation is a better source anyway -- its
# .tbl records the transl_table that was ACTUALLY used, so the re-run reproduces the
# published run rather than re-deriving what it should have been. The class map in
# mito_genetic_codes.json is the fallback, read directly (it is a plain list of
# {code, classes}).
from mito_gene_order import load_origin_anchors, resolve_origin_anchor  # noqa: E402

COX1_PANEL = REPO_ROOT / "assets" / "panels" / "cox1" / "anthozoa.faa"
ANCHORS = REPO_ROOT / "assets" / "taxonomy" / "mito_origin_anchors.json"
GENETIC_CODES = REPO_ROOT / "assets" / "taxonomy" / "mito_genetic_codes.json"

# Only the anthozoan path runs CORAL_ANNOTATION_FIX; every other group grades against its
# reference rather than copying from it, which audit_reference_selection_diff.py covers.
CORAL_CLASSES = {"anthozoa"}

ARMS = ("published", "oldref", "newref")

COLUMNS = ["sample", "old_acc", "old_organism", "new_acc", "new_organism",
           "genetic_code", "origin_gene", "circular", "frame",
           "published_pcgs", "published_nad5_aa", "published_co1_aa",
           "published_internal_stops", "published_gate", "published_sha",
           "oldref_pcgs", "oldref_nad5_aa", "oldref_co1_aa",
           "oldref_internal_stops", "oldref_gate", "oldref_sha",
           "newref_pcgs", "newref_nad5_aa", "newref_co1_aa",
           "newref_internal_stops", "newref_gate", "newref_sha",
           "note", "verdict", "attribution"]


# --------------------------------------------------------------------------- pure logic

def read_fasta(path):
    """{header-id: sequence} from a FASTA, uppercased and whitespace-free."""
    out, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif name is not None:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf).upper()
    return out


def bed_features(path):
    """[(chrom, start, end, name, strand)] from a MITOS result.bed."""
    feats = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track")):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        strand = f[5].strip() if len(f) > 5 else "+"
        feats.append((f[0], int(f[1]), int(f[2]), f[3].strip(), strand))
    return feats


_COMPLEMENT = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")

# MITOS2 gene name -> EMMA gene name, mirroring mitos_to_emma.PCG_MAP and RRNA_MAP.
# Copied rather than imported because mitos_to_emma imports biopython at module level and
# this has to stay importable in a plain python for the unit tests; a test asserts the two
# agree whenever biopython is present, so the copy cannot drift unnoticed.
PCG_MAP = {
    "cox1": "CO1", "cox2": "CO2", "cox3": "CO3", "cob": "CYTB",
    "atp6": "ATP6", "atp8": "ATP8",
    "nad1": "ND1", "nad2": "ND2", "nad3": "ND3",
    "nad4": "ND4", "nad4l": "ND4L", "nad5": "ND5", "nad6": "ND6",
}
RRNA_MAP = {"rrnS": "RNR1", "rrnL": "RNR2"}

# What coral_fix_bed.py rewrites: the 16S rRNA, the intron-split nad5, and cox1 when it is
# intron-split too. These are exactly the genes whose published sequence may no longer
# match the raw MITOS BED, so they cannot anchor a frame check. OG2361 is the case that
# proves it: its raw-BED cox1 is 873 bp while its published CO1 is 1572 bp, because the
# fixer repaired it.
FIXER_TOUCHED = {"CO1", "ND5", "RNR2"}


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def extract(genome_seq, start, end, strand):
    """The nucleotides a BED interval names, oriented to its strand."""
    sub = genome_seq[start:end]
    return revcomp(sub) if strand == "-" else sub


def emma_gene(mitos_name):
    """EMMA gene name for a MITOS BED feature name, or '' if it is not a PCG/rRNA."""
    base = re.split(r"[_(]", (mitos_name or "").strip())[0].lower()
    return PCG_MAP.get(base) or RRNA_MAP.get(base) or RRNA_MAP.get(
        {"rrnl": "rrnL", "rrns": "rrnS"}.get(base, ""), "")


def check_frame(feats, genome_seq, published_cds, min_matches=3):
    """(ok, reason) for a candidate genome against the BED it must match.

    Validated by CONTENT on genes the coral fixer does NOT touch. Two earlier designs
    were wrong against real data and are worth recording:

      * "cox1 sits at BED offset 0" fails the rotate_to_cox1.py fail-safe, which writes an
        assembly through UNROTATED when it finds no confident cox1 hit. MITOS then
        annotates that frame legitimately.
      * "the BED cox1 matches the published CO1" fails whenever the fixer repaired cox1,
        which is precisely the population this audit is about.

    So compare every BED feature that maps to a published CDS, skip the three genes the
    fixer rewrites, and require several independent agreements. A frame that is merely
    shifted matches nothing, so this is a strong test even at min_matches=3.
    """
    if not genome_seq:
        return False, "reconstructed genome is empty"
    if not feats:
        return False, "no BED features"
    longest = max(e for _c, _s, e, _n, _st in feats)
    if longest > len(genome_seq):
        return False, f"BED coordinate {longest} exceeds genome length {len(genome_seq)}"

    matched, mismatched, compared = 0, 0, 0
    for _chrom, start, end, name, strand in feats:
        gene = emma_gene(name)
        if not gene or gene in FIXER_TOUCHED:
            continue
        ref = (published_cds.get(gene) or "").upper()
        if not ref:
            continue
        compared += 1
        got = extract(genome_seq, start, end, strand).upper()
        if got and (got in ref or ref in got):
            matched += 1
        else:
            mismatched += 1
    if compared == 0:
        return False, "no fixer-independent gene available to anchor the frame"
    if matched >= min_matches and matched > mismatched:
        return True, f"ok ({matched}/{compared} untouched genes match)"
    return False, (f"frame does not match the published annotation "
                   f"({matched}/{compared} untouched genes match)")


def published_cds(annot_dir):
    """{EMMA gene: nucleotide sequence} from a published annotation/cds/ directory.

    The filenames already carry the EMMA gene name (CO1.<prefix>.fa), so no mapping is
    needed on this side.
    """
    out = {}
    cds_dir = Path(annot_dir) / "cds"
    if not cds_dir.is_dir():
        return out
    for fa in sorted(cds_dir.glob("*.fa")):
        gene = fa.name.split(".")[0].upper()
        seqs = read_fasta(fa)
        if seqs:
            out[gene] = next(iter(seqs.values()))
    return out


def count_internal_stops(protein_seqs):
    """Stop codons strictly inside a translation, summed over the PCGs.

    A trailing '*' is the stop codon and is correct; anything before it means the reading
    frame is broken, which is the failure a transferred nad5 join can introduce.
    """
    total = 0
    for seq in protein_seqs:
        s = seq.rstrip("*")
        total += s.count("*")
    return total


def arm_metrics(outdir, prefix):
    """PCG count, nad5/CO1 lengths, internal stops and a content signature for one arm.

    The signature is what makes 'identical' meaningful: the re-origined genome plus every
    per-gene translation, so a change in either the deposited sequence or any annotated
    product moves it.
    """
    outdir = Path(outdir)
    prot_dir, cds_dir = outdir / "proteins", outdir / "cds"
    prots = {}
    for fa in sorted(prot_dir.glob("*.fa")) if prot_dir.is_dir() else []:
        gene = fa.name.split(".")[0]
        seqs = read_fasta(fa)
        if seqs:
            prots[gene] = next(iter(seqs.values()))
    genome_fa = outdir / f"{prefix}.fa"
    genome = ""
    if genome_fa.exists():
        seqs = read_fasta(genome_fa)
        genome = next(iter(seqs.values()), "")

    def aa_len(gene_suffix):
        for gene, seq in prots.items():
            if gene.upper().endswith(gene_suffix):
                return len(seq.rstrip("*"))
        return 0

    digest = hashlib.sha256()
    digest.update(genome.encode())
    for gene in sorted(prots):
        digest.update(f"{gene}:{prots[gene]}".encode())

    return {
        "pcgs": len(prots),
        "nad5_aa": aa_len("ND5"),
        "co1_aa": aa_len("CO1"),
        "internal_stops": count_internal_stops(prots.values()),
        "sha": digest.hexdigest()[:16],
        "genome_len": len(genome),
    }


def quality_key(m):
    """Orderable annotation quality: more PCGs and fewer broken frames is better.

    Deliberately coarse. It only has to separate better/worse/neutral, and every component
    is something annotation_qc_gate.py already treats as a defect.
    """
    return (-int(m.get("internal_stops", 0)), int(m.get("pcgs", 0)),
            int(m.get("nad5_aa", 0)), int(m.get("co1_aa", 0)))


def verdict_of(published, oldref, newref):
    """(verdict, attribution) across the three arms.

    Headline is published vs newref, i.e. what changes on disk. Attribution splits that
    difference into the part the code moved and the part the reference moved, which the
    two-arm comparison cannot separate.
    """
    for m in (published, oldref, newref):
        if m is None or m.get("failed"):
            return "failed", "n/a"

    code_moved = published["sha"] != oldref["sha"]
    ref_moved = oldref["sha"] != newref["sha"]
    if code_moved and ref_moved:
        attribution = "both"
    elif ref_moved:
        attribution = "reference"
    elif code_moved:
        attribution = "code_drift"
    else:
        attribution = "n/a"

    if published["sha"] == newref["sha"]:
        return "identical", attribution
    a, b = quality_key(published), quality_key(newref)
    if b > a:
        return "changed_better", attribution
    if b < a:
        return "changed_worse", attribution
    return "changed_neutral", attribution


def load_class_codes(path):
    """{class name (lower): genetic code} from assets/taxonomy/mito_genetic_codes.json."""
    import json
    entries = json.loads(Path(path).read_text()).get("codes", [])
    out = {}
    for entry in entries:
        for name in entry.get("classes", []):
            out[str(name).strip().lower()] = int(entry["code"])
    return out


def published_genetic_code(annot_dir):
    """The transl_table the published annotation actually used, or 0.

    Preferred over re-deriving from class: it is what the run did, so a re-run that
    matches it differs from the published output only where the reference made it differ.
    """
    for tbl in sorted(Path(annot_dir).glob("*.tbl")):
        m = re.search(r"transl_table\s+(\d+)", tbl.read_text())
        if m:
            return int(m.group(1))
    return 0


def circular_from_getorg_check(path):
    """final_verdict_circular from a getorg_check.tsv, defaulting to circular.

    meta.circular == false is what makes the module pass --linear, and a linear genome
    must not be re-origined as if it wrapped. Absent or unparseable, assume circular,
    which is what the pipeline's own default does.
    """
    try:
        rows = list(csv.DictReader(Path(path).open(), delimiter="\t"))
    except OSError:
        return True
    if not rows:
        return True
    value = (rows[0].get("final_verdict_circular") or "").strip().lower()
    if value in ("false", "no", "0"):
        return False
    return True


def reference_accession(gb_path):
    """(accession, organism) from a reference GenBank, without biopython."""
    acc, organism = "", ""
    try:
        text = Path(gb_path).read_text()
    except OSError:
        return acc, organism
    m = re.search(r"^VERSION\s+(\S+)", text, re.M) or re.search(r"^ACCESSION\s+(\S+)", text, re.M)
    if m:
        acc = m.group(1)
    m = re.search(r"^\s*ORGANISM\s+(.+)$", text, re.M)
    if m:
        organism = m.group(1).strip()
    return acc, organism


# ------------------------------------------------------------------------------- runner

def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)


def published_dirs(corpus_dirs):
    """sample -> (prefix, annotation dir, mtdna dir) for every published assembly."""
    found = {}
    for root in corpus_dirs:
        root = Path(root)
        for annot in sorted(root.glob("*/*/annotation")):
            prefix_dir = annot.parent
            sample = prefix_dir.parent.name
            mtdna = prefix_dir / "mtdna"
            if (annot / "mitos_raw" / "result.bed").exists():
                found[sample] = (prefix_dir.name, annot, mtdna)
    return found


def annotate(workdir, tag, raw_bed, rotated, ref_gb, code, prefix, species,
             origin_gene, linear):
    """One arm: patch the BED against ref_gb, then re-run the EMMA adapter."""
    arm_dir = Path(workdir) / tag
    arm_dir.mkdir(parents=True, exist_ok=True)
    fixed_bed = arm_dir / "result.fixed.bed"
    status = arm_dir / f"{prefix}.coral_fix.status.txt"
    fix = run([sys.executable, BIN_DIR / "coral_fix_bed.py",
               "--bed", raw_bed, "--genome", rotated, "--ref-gb", ref_gb,
               "--code", code, "--out-bed", fixed_bed, "--status", status])
    if fix.returncode != 0 or not fixed_bed.exists():
        return None, f"coral_fix_bed failed: {(fix.stderr or '').strip().splitlines()[-1:]}"
    out = arm_dir / "annotation"
    cmd = [sys.executable, BIN_DIR / "mitos_to_emma.py",
           "--bed", fixed_bed, "--genome", rotated, "--prefix", prefix,
           "--outdir", out, "--code", code, "--species", species,
           "--origin-gene", origin_gene]
    if linear:
        cmd.append("--linear")
    emma = run(cmd)
    if emma.returncode != 0:
        return None, f"mitos_to_emma failed: {(emma.stderr or '').strip().splitlines()[-1:]}"
    metrics = arm_metrics(out, prefix)
    metrics["status"] = status.read_text().split("\t")[0].strip() if status.exists() else ""
    metrics["outdir"] = str(out)
    return metrics, ""


def gate(outdir, prefix, code, workdir, tag):
    """annotation_qc_gate.py's verdict for one arm, or '' when it could not run."""
    outdir = Path(outdir)
    gff = outdir / f"{prefix}.gff"
    if not gff.exists():
        return ""
    out = Path(workdir) / f"{tag}.gate.txt"
    res = run([sys.executable, BIN_DIR / "annotation_qc_gate.py",
               "--gff", gff, "--proteins", outdir / "proteins", "--cds", outdir / "cds",
               "--genetic-code", code, "--out", out])
    if res.returncode != 0 or not out.exists():
        return ""
    return out.read_text().split("\t")[0].strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", type=Path, action="append", required=True,
                    help="Pipeline samplesheet CSV. Repeatable.")
    ap.add_argument("--corpus", type=Path, action="append", required=True,
                    help="A published mitogenomes/ directory. Repeatable.")
    ap.add_argument("--refdb-dir", type=Path,
                    default=REPO_ROOT / "assets" / "refdb" / "anthozoa",
                    help="The group database to resolve the NEW reference from.")
    ap.add_argument("--group", default="anthozoa")
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    samples = {}
    for sheet in args.samplesheet:
        for row in csv.DictReader(sheet.open(newline="")):
            if (row.get("sample") or "").strip():
                samples[row["sample"].strip()] = row
    corpus = published_dirs(args.corpus)
    class_codes = load_class_codes(GENETIC_CODES)
    anchors = load_origin_anchors(ANCHORS)
    args.workdir.mkdir(parents=True, exist_ok=True)

    rows, tally, attrib, skipped = [], collections.Counter(), collections.Counter(), collections.Counter()
    for sample in sorted(corpus):
        meta = samples.get(sample)
        if meta is None:
            skipped["no_samplesheet_row"] += 1
            continue
        tax_class = (meta.get("class") or "").strip()
        if tax_class.lower() not in CORAL_CLASSES:
            skipped["not_coral_fix_eligible"] += 1
            continue
        prefix_name, annot, mtdna = corpus[sample]
        raw_bed = annot / "mitos_raw" / "result.bed"
        old_refs = sorted(annot.glob("mitos_fix/*.reference.gb"))
        if not old_refs:
            skipped["no_coral_fix_reference"] += 1   # a PASS coral never ran the fixer
            continue
        assemblies = [p for p in mtdna.glob(f"{prefix_name}.fasta") if p.stat().st_size]
        if not assemblies:
            skipped["no_assembly"] += 1
            continue
        assembly, old_ref = assemblies[0], old_refs[0]
        mitos_prefix = old_ref.name[:-len(".reference.gb")]

        code = published_genetic_code(annot) or class_codes.get(tax_class.lower(), 0)
        if not code:
            skipped["no_genetic_code"] += 1
            continue
        origin_gene = resolve_origin_anchor((meta.get("order") or "").strip(),
                                            tax_class, *anchors)
        checks = sorted(mtdna.glob("*.getorg_check.tsv"))
        circular = circular_from_getorg_check(checks[0]) if checks else True
        work = args.workdir / sample
        work.mkdir(parents=True, exist_ok=True)

        # Reconstruct the frame the BED is in, then prove it before annotating in it.
        # Two candidate frames, verified rather than assumed. ROTATE_ORIGIN has drifted
        # since these corals ran (OG2361's published BED is unrotated, while today's
        # rotate_to_cox1.py does rotate it), so re-running the rotation is not guaranteed
        # to reproduce the frame MITOS actually annotated. Whichever candidate agrees with
        # the published annotation IS the original frame; if neither does, the sample
        # fails rather than being annotated in a frame nobody checked.
        rotated = work / f"{mitos_prefix}.rotated.fa"
        rot = run([sys.executable, BIN_DIR / "rotate_to_cox1.py",
                   "--genome", assembly, "--cox1-ref", COX1_PANEL, "--out", rotated])
        feats = bed_features(raw_bed)
        pub_cds = published_cds(annot)
        candidates = []
        if rot.returncode == 0 and rotated.exists():
            candidates.append(("rotated", rotated))
        candidates.append(("unrotated", assembly))

        note, genome_path, frame_used, frame_reason = "", None, "", ""
        for label, path in candidates:
            seq = next(iter(read_fasta(path).values()), "")
            ok, reason = check_frame(feats, seq, pub_cds)
            if ok:
                genome_path, frame_used, frame_reason = path, label, reason
                break
            frame_reason = reason
        if genome_path is None:
            note = f"frame check failed: {frame_reason}"
        else:
            rotated = genome_path

        new_ref = work / "new.reference.gb"
        if not note:
            sel = run([sys.executable, BIN_DIR / "select_reference_db.py",
                       "--assembly", assembly, "--refdb-dir", args.refdb_dir,
                       "--group", args.group, "--out-gb", new_ref,
                       "--out-status", work / "new.reference_select.txt"])
            if sel.returncode != 0 or not new_ref.exists() or not new_ref.stat().st_size:
                note = "select_reference_db produced no reference"

        species = (meta.get("nominal_species_id") or "").strip()
        arms = {}
        if not note:
            for tag, ref in (("oldref", old_ref), ("newref", new_ref)):
                m, err = annotate(work, tag, raw_bed, rotated, ref, code, mitos_prefix,
                                  species, origin_gene, not circular)
                if m is None:
                    note = err
                    break
                m["gate"] = gate(m["outdir"], mitos_prefix, code, work, tag)
                arms[tag] = m
        if not note:
            arms["published"] = arm_metrics(annot, mitos_prefix)
            arms["published"]["gate"] = gate(annot, mitos_prefix, code, work, "published")

        if note:
            verdict, attribution = "failed", "n/a"
            for tag in ARMS:
                arms.setdefault(tag, {})
        else:
            verdict, attribution = verdict_of(arms["published"], arms["oldref"],
                                              arms["newref"])
        tally[verdict] += 1
        attrib[attribution] += 1
        old_acc, old_org = reference_accession(old_ref)
        new_acc, new_org = reference_accession(new_ref) if new_ref.exists() else ("", "")
        row = {"sample": sample, "old_acc": old_acc, "old_organism": old_org,
               "new_acc": new_acc, "new_organism": new_org,
               "genetic_code": code, "origin_gene": origin_gene,
               "circular": "yes" if circular else "no",
               "frame": frame_used,
               "note": note, "verdict": verdict, "attribution": attribution}
        for tag in ARMS:
            m = arms.get(tag) or {}
            for field in ("pcgs", "nad5_aa", "co1_aa", "internal_stops", "gate", "sha"):
                row[f"{tag}_{field}"] = m.get(field, "")
        rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"[coraldiff] {len(rows)} corals audited -> {args.out}")
    for verdict in ("identical", "changed_neutral", "changed_better", "changed_worse",
                    "failed"):
        if tally[verdict]:
            print(f"[coraldiff]   {verdict:16s} {tally[verdict]}")
    for key, n in sorted(attrib.items()):
        print(f"[coraldiff]   attribution {key:11s} {n}")
    for reason, n in sorted(skipped.items()):
        print(f"[coraldiff]   skipped {reason}: {n}", file=sys.stderr)
    # Name every row that is not a clean 'identical', so a regression cannot hide in a
    # count. changed_worse is the one that would block adopting the rebuild.
    for row in rows:
        if row["verdict"] != "identical":
            print(f"[coraldiff]   {row['verdict']:15s} {row['sample']:10s} "
                  f"{row['old_organism']} -> {row['new_organism']} "
                  f"[{row['attribution']}] {row['note']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
