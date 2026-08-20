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
    p.add_argument("--genetic-code", type=int, default=2,
                   help="NCBI mitochondrial translation table (mgcode). Per-sample value "
                        "from meta.genetic_code: 2=vertebrate, 4=coelenterate/coral, "
                        "9=echinoderm/flatworm (default: 2)")
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
    elif letters.startswith("oatk"):
        asm_prog = "Oatk"
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

def process_fasta_file(input_file, output_file, species, assembly, og_id, genetic_code=2):
    log(f"📝 Processing FASTA file: {input_file}")
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                f_out.write(f">{assembly} [organism={species}] [isolate={og_id}] [mgcode={genetic_code}] {species} mitochondrion\n")
            else:
                f_out.write(line)
    log(f"✅ Processed FASTA: {output_file}")

# Stop codons by NCBI mitochondrial translation table. Used to confirm that an
# incomplete terminal codon really is a truncated stop before we assert aa:TERM.
# Only the vertebrate code (2) reads AGA/AGG as stops; the coelenterate (4),
# invertebrate (5) and echinoderm/flatworm (9) codes use TAA/TAG only (AGA/AGG
# there are Arg or Ser). TGA is Trp, not a stop, in every mitochondrial code.
_MITO_STOPS_BY_CODE = {
    2: ("TAA", "TAG", "AGA", "AGG"),
}
_DEFAULT_MITO_STOPS = ("TAA", "TAG")

def mito_stop_codons(genetic_code: int):
    return _MITO_STOPS_BY_CODE.get(int(genetic_code), _DEFAULT_MITO_STOPS)

# Initiation codons each mitochondrial code accepts without comment. A CDS
# starting on anything else is what raises SEQ_FEAT.StartCodon, so these are the
# codons for which no transl_except is needed. The vertebrate code (2) is the
# strictest and the only one this pipeline routinely uses.
_MITO_STARTS_BY_CODE = {
    2: ("ATT", "ATC", "ATA", "ATG", "GTG"),
}
_DEFAULT_MITO_STARTS = ("ATG", "GTG")

def mito_start_codons(genetic_code: int):
    return _MITO_STARTS_BY_CODE.get(int(genetic_code), _DEFAULT_MITO_STARTS)

# Union of the initiation codons across every mitochondrial translation table,
# read out of the EMBOSS EGC.* data files (tables 2, 3, 4, 5, 9, 13, 14, 21;
# the later 24 and 33 add TTG/CTG/ATG/GTG, all already here).  A CDS starting on
# one of these is initiating on a codon that is a documented start *somewhere*
# in mitochondrial biology, so an alternative-initiation transl_except is a fair
# reading even when the sample's own table does not list it.  A start outside
# this set is far more likely a mis-called CDS boundary or a base-calling error
# than a real alternative start, so those are left to fail table2asn validation
# and get looked at by hand rather than being silently declared as Met.
PLAUSIBLE_MITO_STARTS = ("TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG")
# Substring of EMMA's note ("...TAA stop codon is completed by the addition of
# 3' A residues to the mRNA") that survives the 'putative ' strip. Its presence
# marks a CDS whose stop is completed post-transcriptionally by polyadenylation.
POLYA_NOTE_MARK = "3' A residues"
# Substring of EMMA's note ("non-standard start codon TTG"). Its presence marks a
# CDS that initiates on a codon the genetic code does not list as a start.
START_NOTE_MARK = "non-standard start codon"
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

def compute_transl_except_pos(start: int, end: int, seq: str, stops):
    """For a single-interval CDS (1-based; start>end means minus strand), return
    the transl_except 'pos:' expression if its terminal codon is an incomplete,
    polyadenylation-completed stop, else None. `stops` is the stop-codon tuple for
    this sample's genetic code (see mito_stop_codons).

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
        if not any(stop.startswith(codon) for stop in stops):
            return None
        return f"{hi}" if rem == 1 else f"{lo}..{hi}"
    else:                                    # minus strand; 3' end at low coord
        lo, hi = end, end + rem - 1
        if lo < 1 or hi > n:
            return None
        codon = _revcomp(seq[lo - 1:hi])     # coding-strand partial bases
        if not any(stop.startswith(codon) for stop in stops):
            return None
        return f"complement({lo})" if rem == 1 else f"complement({lo}..{hi})"

def compute_start_transl_except_pos(start: int, end: int, seq: str, starts):
    """Mirror of compute_transl_except_pos for the 5' end. For a single-interval
    CDS (1-based; start>end means minus strand), return (pos, codon) where pos is
    the transl_except 'pos:' expression covering the initiation codon and codon is
    that codon on the coding strand, when the codon is not one the genetic code
    accepts as a start; else None. `starts` is the start-codon tuple for this
    sample's genetic code (see mito_start_codons).

    table2asn raises SEQ_FEAT.StartCodon when a complete CDS begins on a codon
    outside the table's start set, and then emits a gap symbol rather than an
    amino acid for residue 1, which in turn raises SEQ_INST.BadProteinStart.
    Declaring the codon as aa:Met is the INSDC way to record a genuine
    alternative initiation codon, and clears both. Marking the CDS 5'-partial
    would also silence table2asn but would misrepresent a complete gene.

    The codon comes back with the position so the caller can decide whether it is
    a plausible initiator at all (see PLAUSIBLE_MITO_STARTS) and can name it in
    the note it writes.
    """
    n = len(seq)
    if start <= end:                         # plus strand; 5' end at low coord
        lo, hi = start, start + 2
        if lo < 1 or hi > n:
            return None
        codon = seq[lo - 1:hi]               # coding-strand bases
        if codon in starts:
            return None
        return f"{lo}..{hi}", codon
    else:                                    # minus strand; 5' end at high coord
        hi, lo = start, start - 2
        if lo < 1 or hi > n:
            return None
        codon = _revcomp(seq[lo - 1:hi])     # coding-strand bases
        if codon in starts:
            return None
        return f"complement({lo}..{hi})", codon

def product_of(qual_lines):
    """Name a CDS by its /product qualifier, for log messages. Falls back to the
    gene name and then to a placeholder, since neither is guaranteed present."""
    for key in ('product', 'gene'):
        for q in qual_lines:
            cols = q.split('\t')
            if len(cols) >= 2 and cols[-2] == key:
                return cols[-1]
    return 'CDS'

def _tbl_interval(line: str):
    """Return (start, end) if line opens a feature interval, else None.

    Feature lines are non-indented and carry the interval in the first two
    columns.  A feature-key column marks the start of a new feature; a bare
    two-column line is an additional interval of the feature above it.
    """
    if not line or line.startswith('\t'):
        return None
    cols = line.split('\t')
    if len(cols) < 2:
        return None
    try:
        return int(cols[0].lstrip('<>')), int(cols[1].lstrip('<>'))
    except ValueError:
        return None

def sort_tbl_features(lines):
    """Order feature blocks by position on the molecule.

    Emma sorts its output by the start coordinate as a *string*, so a molecule
    comes out as 1, 10053, 1027, 10343, 1099 ...  Re-sort numerically here, at
    the point the feature table is normalised, so the .tbl and everything
    derived from it agree with the coordinate order table2asn and seqret impose
    on the flatfile anyway.  File order is also load bearing downstream: the
    submission pipeline that allocates locus tags numbers loci by walking the
    feature table, so whatever order leaves here is the order the published tags
    carry.

    A block is the feature line, any additional interval lines belonging to a
    joined feature, and the indented qualifier lines that follow.  The sort is
    stable and keyed on the block's lowest coordinate only, which keeps a
    gene adjacent to the mRNA/CDS/tRNA Emma already emits beneath it.
    """
    preamble, blocks = [], []
    for line in lines:
        interval = _tbl_interval(line)
        # A third column holds the feature key and opens a new block; a bare
        # two-column line is a continuation interval of the block above it.
        cols = line.split('\t')
        starts_feature = interval is not None and len(cols) >= 3 and bool(cols[2])
        if starts_feature:
            blocks.append([line])
        elif blocks:
            blocks[-1].append(line)
        else:
            preamble.append(line)
    blocks.sort(key=lambda block: min(_tbl_interval(block[0])))
    return preamble + [line for block in blocks for line in block]

def process_tbl_gb_file(input_file, output_file, assembly, seq=None, genetic_code=2):
    log(f"📝 Processing TBL/GB file: {input_file}")
    stops = mito_stop_codons(genetic_code)
    starts = mito_start_codons(genetic_code)

    def clean(line: str) -> str:
        if line.startswith('>Feature'):
            return f">Feature {assembly}"
        return line.replace('MT-', '').replace('putative ', '')

    with open(input_file) as f_in:
        lines = f_in.read().splitlines()

    out, added, added_start = [], 0, 0
    # CDS whose start EMMA never commented on, and CDS left failing validation
    # because their start codon is not a plausible initiator; both are logged so
    # they are visible without reading the .val.
    unnoted_start, implausible_start = 0, 0
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
            has_start_note = any(START_NOTE_MARK in q for q in qual)
            has_te_met = any('aa:Met' in q for q in qual)
            ic = interval.split('\t')
            # A CDS already marked 5'-partial has no initiation codon to declare;
            # table2asn accepts the truncated start on its own.
            partial_start = ic[0].startswith('<')
            # Only single-interval CDS are handled (mitochondrial genes are single-exon);
            # joined features are passed through untouched to avoid mis-locating the codon.
            # EMMA only sometimes notes a non-standard start, so the note cannot
            # gate this -- the sequence itself decides. What the note does not
            # tell us either way is whether the codon is a credible initiator, so
            # that is checked against PLAUSIBLE_MITO_STARTS instead.
            if (seq is not None and not has_te_met
                    and not partial_start and not extra_intervals):
                found = compute_start_transl_except_pos(int(ic[0].lstrip('<>')),
                                                        int(ic[1].lstrip('<>')), seq, starts)
                if found:
                    pos, codon = found
                    if codon in PLAUSIBLE_MITO_STARTS:
                        out.append(f"\t\t\ttransl_except\t(pos:{pos},aa:Met)")
                        added_start += 1
                        if not has_start_note:
                            out.append(f"\t\t\tnote\t{START_NOTE_MARK} {codon}")
                            unnoted_start += 1
                    else:
                        implausible_start += 1
                        log(f"⚠️  {product_of(qual)} starts on {codon}, which no "
                            f"mitochondrial code lists as an initiator -- leaving it "
                            f"for table2asn to flag rather than declaring aa:Met")
            if seq is not None and has_note and not has_te and not extra_intervals:
                pos = compute_transl_except_pos(int(ic[0].lstrip('<>')),
                                                int(ic[1].lstrip('<>')), seq, stops)
                if pos:
                    out.append(f"\t\t\ttransl_except\t(pos:{pos},aa:TERM)")
                    added += 1
            for q in qual:
                # EMMA's protein_id is an internal placeholder UUID
                # (gnl|Emma|<uuid>), not a real INSDC accession. ENA/GenBank
                # assign the real protein_id at accessioning time, so a
                # submitter-supplied value here is invalid and gets rejected
                # (EMBOSS demotes it to a stray /note on EMBL conversion).
                if '\tprotein_id\t' in q:
                    continue
                out.append(clean(q))
            i = j
            continue
        out.append(clean(lines[i]))
        i += 1

    out = sort_tbl_features(out)

    with open(output_file, 'w') as f_out:
        f_out.write('\n'.join(out) + '\n')
    log(f"✅ Processed TBL/GB: {output_file} "
        f"(transl_except aa:TERM added to {added} CDS, aa:Met to {added_start} CDS, "
        f"{unnoted_start} of those unnoted by EMMA; "
        f"{implausible_start} CDS left with an implausible start codon)")

def sort_gff_records(lines):
    """Order GFF loci by position, mirroring sort_tbl_features.

    Emma's string sort affects the GFF the same way it affects the .tbl.  A
    block is a top-level record plus the records that hang off it via
    Parent=, so gene -> mRNA -> CDS stays together even where a child's
    coordinates differ from its parent's.  Directives and the whole-molecule
    'region' record are pinned ahead of the sorted blocks; region spans the
    entire sequence, so sorting it by start would let the first gene overtake
    it.  Unlike the .tbl, GFF always writes start <= end with the strand in
    column 7, so column 4 alone is the key.
    """
    preamble, blocks = [], []
    for line in lines:
        fields = line.rstrip('\r\n').split('\t')
        record = len(fields) == 9 and not line.startswith('#')
        if record and fields[2] != 'region' and 'Parent=' not in fields[8]:
            blocks.append([line])
        elif blocks:
            blocks[-1].append(line)
        else:
            preamble.append(line)
    blocks.sort(key=lambda block: int(block[0].split('\t')[3]))
    return preamble + [line for block in blocks for line in block]

def process_gff_file(input_file, output_file, assembly):
    log(f"📝 Processing GFF file: {input_file}")
    with open(input_file, 'r') as f:
        original_lines = f.readlines()

    new_lines = []

    for line in original_lines:
        if line.startswith('##sequence-region'):
            cols = line.rstrip('\r\n').split('\t')
            if len(cols) >= 2:
                cols[1] = assembly
            line = '\t'.join(cols) + '\n'
            new_lines.append(line)
            continue

        if line.startswith('#') or not line.strip():
            new_lines.append(line)
            continue

        fields = line.rstrip('\r\n').split('\t')
        if len(fields) == 9:
            fields[0] = assembly
            attrs = fields[8]
            for old, new in NAME_MAP.items():
                attrs = attrs.replace(old, new)
            attrs = re.sub(r';{2,}', ';', attrs).strip(';')
            fields[8] = attrs

        processed_line = '\t'.join(fields)
        processed_line = re.sub(r'(?i)\bputative\b[\s;,:]*', '', processed_line)
        processed_line = processed_line.replace('MT-', '')
        new_lines.append(processed_line + '\n')

    new_lines = sort_gff_records(new_lines)

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

    GENETIC_CODE = args.genetic_code
    process_fasta_file(fasta_in, fasta_out, SPECIES, ASSEMBLY, OG_ID, genetic_code=GENETIC_CODE)
    tmp_gff = outdir / (gff_out.name + ".tmp")
    process_gff_file(gff_in, tmp_gff, ASSEMBLY)
    write_meta_directives(tmp_gff, SPECIES, OG_ID)
    tmp_gff.replace(gff_out)
    # Sequence is needed to resolve polyA-completed stop codons into transl_except
    # qualifiers (clears table2asn SEQ_FEAT.NoStop). Header rewriting does not touch
    # the bases, so coordinates match either the input or output FASTA.
    seq = read_fasta_sequence(fasta_in)
    process_tbl_gb_file(gb_in, gb_out, ASSEMBLY, seq=seq, genetic_code=GENETIC_CODE)

    log("🏁 File processing complete.")

if __name__ == "__main__":
    main()
