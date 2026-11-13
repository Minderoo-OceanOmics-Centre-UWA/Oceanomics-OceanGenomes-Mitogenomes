#!/usr/bin/env python3
"""Extend EMMA protein-coding genes whose stop codon EMMA cut short.

EMMA sometimes ends a PCG one codon early and annotates the last SENSE codon as
a "putative ... stop codon completed by the addition of 3' A residues", even
though a complete TAA/TAG is encoded immediately downstream. Whether it does so
depends on where the assembler linearised the circle, not on the sequence:
NOVA_260925_LA had four Glaucosoma hebraicum assemblies that are the SAME
genome as a passing batch-mate, only rotated (OG3104/OG3157 vs OG3208 at
+12,319 bp; OG3193/OG3250 vs their twins at +9,087 bp). EMMA's --rotate MT-TF
then hides the difference in the published files. Nine had CO1 ending
"...GCA AAA|TAA" and two had ND1 ending "...CTC|TAA"; the CDS integrity check in
annotation_stats.py held all eleven, correctly, and vertebrates had no repair
path (orf_snap / INVERT_ANNOTATION_FIX only see MITOS2 output).

The repair is reference-free and only ever moves a CDS's 3' end downstream by
1-3 nt, onto a full stop the assembly itself encodes:

    span length % 3 == 0, last codon not a stop  -> extend 3 nt
    span ends in a partial stop (T / TA)         -> extend 2 / 1 nt

and only when ALL of these hold:

    * the extended CDS ends on TAA or TAG (AGA/AGG are not accepted as a NEW
      stop: in vertebrate mitogenomes those are usually -1 frameshift sites);
    * the added bases lie inside no other annotated feature, on either strand.
      This is what keeps EMMA's genuine poly-A stops (e.g. COX2 ending in T with
      tRNA-Lys starting right after) exactly as they are;
    * the extension does not cross the linearisation point;
    * the extended CDS has no internal stop, and its translation agrees with
      EMMA's protein over the existing length.

A gene that fails any guard is left as EMMA wrote it, so the sample is held at
the QC gate exactly as before. Edits touch the .gff (gene/mRNA/CDS end, Note
dropped), the .tbl (coordinates, note dropped) and that gene's cds/ and
proteins/ FASTAs. The .svg is not redrawn; it is a picture, not a record.

A status file lists one line per repaired or skipped gene (REPAIRED / SKIP) or
a single NONE line, and the script always exits 0.

Usage:
    emma_stop_repair.py --annotation-dir annotation --code 2 \\
        --status <prefix>.stop_repair.status.txt [--dry-run]
"""

import argparse
import re
import sys
from pathlib import Path

from orf_utils import revcomp, stop_codons, translate

# A complete stop the repair may move a CDS onto. Kept to the two codons that are
# stops in every vertebrate/invertebrate mitochondrial table this pipeline uses.
EXTEND_STOPS = ("TAA", "TAG")

POLY_A_NOTE_RE = re.compile(
    r"putative \S+ stop codon is completed by the addition of 3' A residues to the mRNA")

FASTA_WIDTH = 70


# ------------------------------------------------------------------- parsing

def read_fasta(path):
    header, seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:]
        elif line.strip():
            seq.append(line.strip())
    return header, "".join(seq).upper()


def write_fasta(path, header, seq, width=FASTA_WIDTH):
    with open(path, "w") as out:
        out.write(f">{header}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i:i + width] + "\n")


def parse_attrs(col):
    out = {}
    for part in col.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def read_gff(path):
    """Return (lines, rows). rows: dicts for each 9-column feature line."""
    lines = open(path).read().splitlines()
    rows = []
    for i, line in enumerate(lines):
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        rows.append({"i": i, "type": f[2], "start": int(f[3]), "end": int(f[4]),
                     "strand": f[6], "attrs": parse_attrs(f[8])})
    return lines, rows


def find_one(directory, pattern):
    hits = sorted(Path(directory).glob(pattern))
    return hits[0] if len(hits) == 1 else None


# ------------------------------------------------------------------ planning

def plan_gene(name, cds_row, gene_names_rows, genome, code, prot_seq):
    """Decide whether one CDS needs and permits a stop extension.

    Returns (None, None) when the CDS already ends on a full stop or on a genuine
    poly-A-completed one, (None, reason) when it is short but a guard refused the
    edit, or (edit, None).
    """
    edit, reason = _plan_gene(name, cds_row, gene_names_rows, genome, code, prot_seq)
    if reason and _is_poly_a_partial(cds_row, genome, code):
        # A T / TA tail the assembly does not complete (or completes only inside
        # the next feature) is EMMA's legitimate poly-A convention, not a defect.
        return None, None
    return edit, reason


def _is_poly_a_partial(cds_row, genome, code):
    s, e, strand = cds_row["start"], cds_row["end"], cds_row["strand"]
    span = genome[s - 1:e] if strand == "+" else revcomp(genome[s - 1:e])
    rem = len(span) % 3
    return bool(rem) and any(st.startswith(span[-rem:]) for st in stop_codons(code))


def _plan_gene(name, cds_row, gene_names_rows, genome, code, prot_seq):
    s, e, strand = cds_row["start"], cds_row["end"], cds_row["strand"]
    L = len(genome)
    if s > e or strand not in "+-":
        return None, "unsupported_coordinates"
    span = genome[s - 1:e] if strand == "+" else revcomp(genome[s - 1:e])
    stops = stop_codons(code)
    rem = len(span) % 3
    if rem == 0 and span[-3:] in stops:
        return None, None                      # already a complete stop
    k = 3 - rem                                # 3 for rem 0, else complete the codon

    if strand == "+":
        if e + k > L:
            return None, "extension_crosses_origin"
        added = range(e + 1, e + k + 1)
        new_s, new_e = s, e + k
        new_span = span + genome[e:e + k]
    else:
        if s - k < 1:
            return None, "extension_crosses_origin"
        added = range(s - k, s)
        new_s, new_e = s - k, e
        new_span = revcomp(genome[s - k - 1:e])

    stop = new_span[-3:]
    if stop not in EXTEND_STOPS or stop not in stops:
        return None, f"no_encoded_stop({stop})"

    for other in gene_names_rows:
        if other["name"] == name:
            continue
        if any(other["start"] <= p <= other["end"] for p in added):
            return None, f"overlaps_{other['name']}"

    aa = translate(new_span, code)
    if "*" in aa[:-1]:
        return None, "internal_stop"
    new_cds = new_span[:-3]
    new_aa = translate(new_cds, code)
    if prot_seq is not None:
        # EMMA writes Met for an alternative initiator, so compare from residue 2.
        old = prot_seq.rstrip("*")
        if not old or new_aa[1:len(old)] != old[1:]:
            return None, "protein_mismatch"
        new_prot = old[0] + new_aa[1:]
    else:
        new_prot = None

    return {
        "name": name, "strand": strand,
        "old": (s, e), "new": (new_s, new_e), "k": k, "stop": stop,
        "new_cds": new_cds, "new_prot": new_prot,
    }, None


# ------------------------------------------------------------------- writing

def strip_poly_a_note(attr_col):
    parts = [p for p in attr_col.split(";")
             if not (p.startswith("Note=") and POLY_A_NOTE_RE.search(p))]
    return ";".join(parts)


def apply_gff(lines, rows, edit):
    s, e = edit["old"]
    ns, ne = edit["new"]
    gene = edit["name"]
    for r in rows:
        if r["attrs"].get("Name") != gene or r["type"] not in ("gene", "mRNA", "CDS"):
            continue
        if (r["start"], r["end"]) != (s, e):
            continue
        f = lines[r["i"]].split("\t")
        f[3], f[4] = str(ns), str(ne)
        f[8] = strip_poly_a_note(f[8])
        lines[r["i"]] = "\t".join(f)


def apply_tbl(tbl_lines, edit):
    s, e = edit["old"]
    ns, ne = edit["new"]
    old = (s, e) if edit["strand"] == "+" else (e, s)
    new = (ns, ne) if edit["strand"] == "+" else (ne, ns)
    out, in_cds = [], False
    for line in tbl_lines:
        f = line.split("\t")
        if line and not line.startswith("\t") and len(f) >= 3 and f[0].isdigit():
            in_cds = False
            if (int(f[0]), int(f[1])) == old and f[2] in ("gene", "mRNA", "CDS"):
                f[0], f[1] = str(new[0]), str(new[1])
                line = "\t".join(f)
                in_cds = f[2] == "CDS"
        elif in_cds and line.strip().startswith("note") and POLY_A_NOTE_RE.search(line):
            continue
        out.append(line)
    return out


# -------------------------------------------------------------------- driver

def repair(ann_dir, code, dry_run=False):
    ann = Path(ann_dir)
    gff = find_one(ann, "*.gff")
    tbl = find_one(ann, "*.tbl")
    fa = find_one(ann, "*.fa")
    if not (gff and tbl and fa):
        return [("SKIP", "-", "bundle_incomplete(gff/tbl/fa)")]

    _, genome = read_fasta(fa)
    lines, rows = read_gff(gff)
    gene_rows = [{"name": r["attrs"].get("Name", ""), "start": r["start"], "end": r["end"]}
                 for r in rows if r["type"] == "gene"]

    edits, status = [], []
    for r in rows:
        if r["type"] != "CDS":
            continue
        name = r["attrs"].get("Name", "")
        bare = name[3:] if name.startswith("MT-") else name
        prot_path = find_one(ann / "proteins", f"{name}.*.fa") if (ann / "proteins").is_dir() else None
        prot = read_fasta(prot_path)[1] if prot_path else None
        edit, reason = plan_gene(name, r, gene_rows, genome, code, prot)
        if edit:
            edit["bare"] = bare
            edits.append(edit)
        elif reason:
            status.append(("SKIP", bare, reason))

    for ed in edits:
        (s, e), (ns, ne) = ed["old"], ed["new"]
        status.append(("REPAIRED", ed["bare"],
                       f"stop_extended(+{ed['k']},{ed['stop']}) {s}-{e} -> {ns}-{ne} {ed['strand']}"))

    if dry_run or not edits:
        return status or [("NONE", "-", "no truncated stops")]

    # Build every rewritten text before touching disk, so a failure part-way can
    # never leave a .gff and .tbl that disagree.
    for ed in edits:
        apply_gff(lines, rows, ed)
    tbl_lines = tbl.read_text().splitlines()
    for ed in edits:
        tbl_lines = apply_tbl(tbl_lines, ed)
    gff.write_text("\n".join(lines) + "\n")
    tbl.write_text("\n".join(tbl_lines) + "\n")

    for ed in edits:
        cds_path = find_one(ann / "cds", f"{ed['name']}.*.fa")
        if cds_path:
            write_fasta(cds_path, read_fasta(cds_path)[0], ed["new_cds"])
        prot_path = find_one(ann / "proteins", f"{ed['name']}.*.fa")
        if prot_path and ed["new_prot"]:
            write_fasta(prot_path, read_fasta(prot_path)[0], ed["new_prot"])
    return status


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--code", required=True, type=int)
    ap.add_argument("--status", help="status TSV (default: stdout)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change without writing")
    a = ap.parse_args()

    try:
        status = repair(a.annotation_dir, a.code, a.dry_run)
    except Exception as exc:          # fail-safe: never break the bundle
        status = [("SKIP", "-", f"crashed: {type(exc).__name__}: {exc}")]
    text = "".join("\t".join(row) + "\n" for row in status)
    if a.status:
        Path(a.status).write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
