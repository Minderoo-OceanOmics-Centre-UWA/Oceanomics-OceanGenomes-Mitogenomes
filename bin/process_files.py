#!/usr/bin/env python3
import os
import re
import argparse
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Process mitogenome files for GenBank submission.")
    p.add_argument("--og-id", required=True, help="OG_ID / isolate")
    p.add_argument("--species", required=True, help="Species name")
    # directory in / out
    p.add_argument("--input-dir", required=True, help="Directory containing .fa/.fasta, .gff, and .gb/.tbl files")
    p.add_argument("--outdir", default="processed", help="Output directory (default: processed)")
    return p.parse_args()

def log(msg):
    print(msg, flush=True)

NAME_MAP = {
    "Name=12srna":"Name=RNR1","Name=16srna":"Name=RNR2","Name=trnH-GUG":"Name=TH","Name=trnW-UCA":"Name=TW",
    "Name=trnA-UGC":"Name=TA","Name=trnN-GUU":"Name=TN","Name=trnC-GCA":"Name=TC","Name=trnY-GUA":"Name=TY",
    "Name=trnS1-GCU":"Name=TS1","Name=trnL1-UAG":"Name=TL1","Name=trnS2-UGA":"Name=TS2","Name=trnD-GUC":"Name=TD",
    "Name=trnK-UUU":"Name=TK","Name=trnG-UCC":"Name=TG","Name=trnR-UCG":"Name=TR","Name=trnE-UUC":"Name=TE",
    "Name=trnT-UGU":"Name=TT","Name=trnP-UGG":"Name=TP","Name=trnF-GAA":"Name=TF","Name=trnV-UAC":"Name=TV",
    "Name=trnL2-UAA":"Name=TL2","Name=trnI-GAU":"Name=TI","Name=trnQ-UUG":"Name=TQ","Name=trnM-CAU":"Name=TM",
}

# --- NEW: choose newest matching files in input_dir ---
def newest_match(dirpath: Path, patterns):
    matches = []
    for pat in patterns:
        matches.extend(dirpath.glob(pat))
    return max(matches, key=lambda p: p.stat().st_mtime) if matches else None

def _find_matching(d: Path, stem: str, patterns):
    # Return the first file whose basename (without extension) == stem
    for pat in patterns:
        for p in sorted(d.glob(pat)):
            if p.stem == stem:
                return p
    return None

def resolve_three_inputs(input_dir: str):
    """
    Scan the directory, find all stems that have FASTA + GFF + (TBL or GB),
    and return the triplet for the newest FASTA.
    """
    d = Path(input_dir)
    if not d.is_dir():
        raise SystemExit(f"Input dir not found: {d}")

    fastas = list(d.glob("*.fa")) + list(d.glob("*.fasta"))
    if not fastas:
        raise SystemExit(f"Missing in {d}: FASTA")

    # Build candidate triplets by shared stem (basename without extension).
    candidates = []
    for fa in fastas:
        stem = fa.stem
        gff = _find_matching(d, stem, ["*.gff"])
        feat = _find_matching(d, stem, ["*.tbl", "*.gb"])  # prefer .tbl if present
        if gff and feat:
            candidates.append((fa.stat().st_mtime, fa, gff, feat, stem))

    if not candidates:
        raise SystemExit(
            f"No complete triplet found in {d}. "
            "Need FASTA + GFF + (TBL or GB) sharing the same basename."
        )

    # Choose the triplet whose FASTA is newest
    _, fasta, gff, feat, stem = max(candidates, key=lambda t: t[0])
    return fasta, gff, feat, stem

# --- NEW: replicate your bash cmt derivation ---
def derive_seqid_from_fasta(fa_path: Path) -> str:
    return fa_path.stem  # basename without extension

def derive_assembly_method(seqid: str) -> str:
    parts = seqid.split(".")
    if len(parts) < 4:
        raise SystemExit(f"❌  SeqID '{seqid}' has fewer than 4 dot-fields; cannot derive assembly method")
    field4 = parts[3]                       # e.g. v177getorg
    clean  = field4.lstrip("v")             # -> 177getorg
    digits = "".join(ch for ch in clean if ch.isdigit())  # -> 177
    letters = "".join(ch for ch in clean if ch.isalpha()).lower()  # -> getorg

    version = ".".join(list(digits)) if digits else ""
    if letters.startswith("getorg"):
        asm_prog = "GetOrganelle"
    elif letters.startswith("mitohifi"):
        asm_prog = "MitoHifi"
    else:
        raise SystemExit(f"❌  Unknown assembler code in '{seqid}' (field '{field4}')")
    return f"{asm_prog} v.{version}" if version else asm_prog

def derive_seq_tech(seqid: str) -> str:
    s = seqid.lower()
    if "hifi" in s:
        return "PacBio HiFi"
    if "ilmn" in s:
        return "Illumina"
    if "hic" in s:
        return "Hi-C"
    raise SystemExit(f"❌  Unknown sequencing tech in '{seqid}'")

def write_cmt(outdir: Path, seqid: str, assembly_method: str, seq_tech: str) -> Path:
    cmt_path = outdir / f"{seqid}.cmt"
    with open(cmt_path, "w") as f:
        f.write("SeqID\tStructuredCommentPrefix\tAssembly Method\tSequencing Technology\tStructuredCommentSuffix\n")
        f.write(f"{seqid}\tAssembly-Data\t{assembly_method}\t{seq_tech}\tAssembly-Data\n")
    log(f"✅  Wrote {cmt_path.name}   (assembly = '{assembly_method}' ; tech = '{seq_tech}')")
    return cmt_path

def process_fasta_file(input_file, output_file, species, assembly, og_id):
    log(f"📝 Processing FASTA file: {input_file}")
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(f">{assembly} [organism={species}] [isolate={og_id}] [mgcode=2] {species} mitochondrion\n")
            else:
                f_out.write(line)
    log(f"✅ Processed FASTA: {output_file}")

# Vertebrate mitochondrial stop codons (transl_table 2). Used to confirm that an
# incomplete terminal codon really is a truncated stop before we assert aa:TERM.
VERT_MITO_STOPS = ("TAA", "TAG", "AGA", "AGG")
# Substring of EMMA's note ("...TAA stop codon is completed by the addition of
# 3' A residues to the mRNA") that survives the 'putative ' strip. Its presence
# marks a CDS whose stop is completed post-transcriptionally by polyadenylation.
POLYA_NOTE_MARK = "3' A residues"
_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def read_fasta_sequence(fa_path) -> str:
    seq = []
    with open(fa_path) as f:
        for line in f:
            if not line.startswith('>'):
                seq.append(line.strip())
    return ''.join(seq).upper()

def _revcomp(s: str) -> str:
    return ''.join(_COMP.get(b, 'N') for b in reversed(s))

def compute_transl_except_pos(start: int, end: int, seq: str):
    """For a single-interval CDS (1-based; start>end means minus strand), return
    the transl_except 'pos:' expression if its terminal codon is an incomplete,
    polyadenylation-completed stop, else None.

    table2asn raises SEQ_FEAT.NoStop when a CDS lacks an in-frame stop and cannot
    extend into the genome to find one (because the stop is really created by the
    3' poly(A) tail, not encoded). Declaring the truncated codon as aa:TERM tells
    table2asn the terminator is present, clearing the error.
    """
    length = abs(end - start) + 1
    rem = length % 3
    if rem == 0:
        return None
    n = len(seq)
    if start <= end:                         # plus strand; 3' end at high coord
        hi, lo = end, end - rem + 1
        if lo < 1 or hi > n:
            return None
        codon = seq[lo - 1:hi]               # coding-strand partial bases
        if not any(stop.startswith(codon) for stop in VERT_MITO_STOPS):
            return None
        return f"{hi}" if rem == 1 else f"{lo}..{hi}"
    else:                                    # minus strand; 3' end at low coord
        lo, hi = end, end + rem - 1
        if lo < 1 or hi > n:
            return None
        codon = _revcomp(seq[lo - 1:hi])     # coding-strand partial bases
        if not any(stop.startswith(codon) for stop in VERT_MITO_STOPS):
            return None
        return f"complement({lo})" if rem == 1 else f"complement({lo}..{hi})"

def process_tbl_gb_file(input_file, output_file, assembly, seq=None):
    log(f"📝 Processing TBL/GB file: {input_file}")

    def clean(line: str) -> str:
        if line.startswith('>Feature'):
            return f">Feature {assembly}"
        return line.replace('MT-', '').replace('putative ', '')

    with open(input_file) as f_in:
        lines = f_in.read().splitlines()

    out, added = [], 0
    i, n = 0, len(lines)
    while i < n:
        cols = lines[i].split('\t')
        # A CDS feature starts on a non-indented line: <start> <end> CDS
        if len(cols) >= 3 and cols[2] == 'CDS' and not lines[i].startswith('\t'):
            interval = lines[i]
            j = i + 1
            # Additional interval lines (joined CDS) are non-indented "<start> <end>".
            extra_intervals = []
            while (j < n and not lines[j].startswith('\t')
                   and len(lines[j].split('\t')) == 2
                   and lines[j].split('\t')[0].lstrip('<>').isdigit()):
                extra_intervals.append(lines[j]); j += 1
            # Qualifier lines are indented.
            qual = []
            while j < n and lines[j].startswith('\t'):
                qual.append(lines[j]); j += 1

            out.append(clean(interval))
            for e in extra_intervals:
                out.append(clean(e))

            has_note = any(POLYA_NOTE_MARK in q for q in qual)
            has_te = any('\ttransl_except\t' in q or q.lstrip().startswith('transl_except') for q in qual)
            # Only single-interval CDS are handled (mitochondrial genes are single-exon);
            # joined features are passed through untouched to avoid mis-locating the stop.
            if seq is not None and has_note and not has_te and not extra_intervals:
                ic = interval.split('\t')
                pos = compute_transl_except_pos(int(ic[0].lstrip('<>')),
                                                int(ic[1].lstrip('<>')), seq)
                if pos:
                    out.append(f"\t\t\ttransl_except\t(pos:{pos},aa:TERM)")
                    added += 1
            for q in qual:
                out.append(clean(q))
            i = j
            continue
        out.append(clean(lines[i]))
        i += 1

    with open(output_file, 'w') as f_out:
        f_out.write('\n'.join(out) + '\n')
    log(f"✅ Processed TBL/GB: {output_file} (transl_except added to {added} CDS)")

def process_gff_file(input_file, output_file, assembly):
    log(f"📝 Processing GFF file: {input_file}")
    with open(input_file, 'r') as f:
        original_lines = f.readlines()

    new_lines = []
    in_header_block = True

    for line in original_lines:
        if in_header_block and line.startswith('##'):
            if line.startswith('##sequence-region'):
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 2:
                    cols[1] = assembly
                line = '\t'.join(cols) + '\n'
            new_lines.append(line)
            continue

        if in_header_block:
            new_lines.append(line)
            in_header_block = False
            continue

        if line.startswith('#'):
            new_lines.append(line)
            continue

        fields = line.rstrip('\n').split('\t')
        if fields:
            fields[0] = assembly
            if len(fields) > 8:
                attrs = fields[8]
                for old, new in NAME_MAP.items():
                    attrs = attrs.replace(old, new)
                attrs = re.sub(r';{2,}', ';', attrs).strip(';')
                fields[8] = attrs

            processed_line = '\t'.join(fields)
            processed_line = re.sub(r'(?i)\bputative\b[\s;,:]*', '', processed_line)
            processed_line = processed_line.replace('MT-', '')
            new_lines.append(processed_line + '\n')

    with open(output_file, 'w') as f:
        for line in new_lines:
            clean = re.sub(r'(?i)\bputative\b[\s;,:]*', '', line)
            clean = clean.replace('MT-', '')
            f.write(clean)
    log(f"✅ Processed GFF: {output_file}")

def write_meta_directives(output_gff, species, og_id):
    with open(output_gff, 'r') as f:
        lines = f.readlines()
    result = []
    i = 0
    while i < len(lines) and lines[i].startswith('##'):
        result.append(lines[i]); i += 1
    result.append(f"##organism {species}\n")
    result.append(f"##isolate {og_id}\n")
    result.extend(lines[i:])
    with open(output_gff, 'w') as f:
        f.writelines(result)

def main():
    args = parse_args()
    OG_ID = args.og_id
    SPECIES = args.species


    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # find exactly one fasta, gff, gb/tbl in the input directory
    fasta_in, gff_in, gb_in, matched_stem = resolve_three_inputs(args.input_dir)
    ASSEMBLY = matched_stem  # override CLI to mirror your staging planner logic

    # --- Structured Comment: use the SAME SeqID as the FASTA header (== ASSEMBLY) ---
    assembly_method = derive_assembly_method(ASSEMBLY)
    seq_tech = derive_seq_tech(ASSEMBLY)
    write_cmt(outdir, ASSEMBLY, assembly_method, seq_tech)

    log(f"🔍 Sample: {OG_ID}\n   Species: {SPECIES}\n   Assembly: {ASSEMBLY}")
    log(f"📁 Inputs:\n   FASTA: {fasta_in}\n   GFF:   {gff_in}\n   GB/TBL:{gb_in}")
    log(f"📦 Output dir: {outdir}")

    # Output paths keep original basenames
    fasta_out = outdir / Path(fasta_in).name
    gff_out   = outdir / Path(gff_in).name
    gb_out    = outdir / Path(gb_in).name

    process_fasta_file(fasta_in, fasta_out, SPECIES, ASSEMBLY, OG_ID)
    tmp_gff = outdir / (gff_out.name + ".tmp")
    process_gff_file(gff_in, tmp_gff, ASSEMBLY)
    write_meta_directives(tmp_gff, SPECIES, OG_ID)
    tmp_gff.replace(gff_out)
    # Sequence is needed to resolve polyA-completed stop codons into transl_except
    # qualifiers (clears table2asn SEQ_FEAT.NoStop). Header rewriting does not touch
    # the bases, so coordinates match either the input or output FASTA.
    seq = read_fasta_sequence(fasta_in)
    process_tbl_gb_file(gb_in, gb_out, ASSEMBLY, seq=seq)

    log("🏁 File processing complete.")

if __name__ == "__main__":
    main()
