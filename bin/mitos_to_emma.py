#!/usr/bin/env python3
"""
Adapt MITOS2 annotation output to the EMMA output contract used by the rest of
the OceanGenomes mitogenome pipeline.

MITOS2 emits result.bed / result.gff / result.faa / result.fas with its own gene
nomenclature (cox1, cob, rrnS, rrnL, trnF(gaa), ...). Downstream modules
(annotation_stats.py, extract_genes_gff.py, process_files.py, the BLAST->LCA mix)
were written against EMMA, which uses MT-<GENE> naming (CO1, CYTB, RNR1, RNR2,
TF, ...), an EMMA-style GFF (##sequence-region + gene/CDS/tRNA/rRNA lines with
Name= attributes) and per-gene FASTAs under cds/ and proteins/.

This script reads MITOS2's result.bed (stable, well-defined columns) plus the
assembly FASTA and writes an EMMA-shaped directory:

    <outdir>/<prefix>.fa                      (the assembly genome)
    <outdir>/<prefix>.gff                     (EMMA-style GFF)
    <outdir>/cds/<GENE>.<prefix>.fa           (per-gene nucleotide FASTA)
    <outdir>/proteins/MT-<GENE>.<prefix>.fa   (per-PCG translated FASTA)

Coordinates and sequences are derived from the BED so the result is independent
of MITOS2's GFF/FASTA header formatting, which varies between releases.

Re-origin to --origin-gene: the EMMA-shaped genome is published starting at the
5' end of the anchor gene the caller names (reverse-complementing first if that
gene is on the minus strand). The anchor's position comes from MITOS's own
annotation -- no extra tool or BLAST is needed -- and the genome, GFF and per-gene
FASTAs are all shifted together so they stay consistent.

The anchor is per-taxon, resolved by the caller from
assets/taxonomy/mito_origin_anchors.json. It is NOT trnM for everything, which is
what this script used to do while claiming trnM matched "how coral mitogenomes are
deposited in NCBI". Measured across the curated RefSeq databases in assets/refdb/,
trnM is the deposited origin for Porifera 0%, Annelida 0%, Ctenophora 0%,
Echinodermata 0.7%, Mollusca 0.8% and Arthropoda 1.9%. It IS the convention for
Scleractinia (63.8%), which is where the rule came from and where it still applies;
Anthozoa as a whole is rrnL 33.9% / cox1 28.5% / trnM 17.6%, i.e. no majority.

Two rotations, two different jobs. The upstream cox1 pre-rotation
(rotate_to_cox1.py) runs first purely to move the LINEARISATION POINT off an
intron-split gene so MITOS annotates cleanly, and cox1 is a good universal choice
for that because it is conserved and findable by tblastn. This step then sets the
PUBLISHED origin, which is a deposition convention and has to be per-taxon. The
second cannot be folded into the first: tRNA-Met and rrnL are not findable by
protein search, which is why they need MITOS's annotation to exist first.

Skipped for linear molecules. When the anchor gene is not annotated the genome
falls back to --origin-fallback (cox1 by default, which is where the pre-rotation
already put it) and, failing that, stays on the pre-rotation origin.

Usage:
    mitos_to_emma.py --bed result.bed --genome asm.fa \
        --prefix OG1.ilmn.250131.getorg1770.mitos2110 --outdir emma --code 5
"""

import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# MITOS2 protein-coding gene name -> EMMA gene name
PCG_MAP = {
    "cox1": "CO1", "cox2": "CO2", "cox3": "CO3",
    "cob": "CYTB",
    "atp6": "ATP6", "atp8": "ATP8",
    "nad1": "ND1", "nad2": "ND2", "nad3": "ND3",
    "nad4": "ND4", "nad4l": "ND4L", "nad5": "ND5", "nad6": "ND6",
}

# MITOS2 rRNA name -> EMMA gene name
RRNA_MAP = {"rrnS": "RNR1", "rrnL": "RNR2"}

# Products for nicer GenBank output (downstream falls back to gene name otherwise)
PRODUCT = {
    "CO1": "cytochrome c oxidase subunit I",
    "CO2": "cytochrome c oxidase subunit II",
    "CO3": "cytochrome c oxidase subunit III",
    "CYTB": "cytochrome b",
    "ATP6": "ATP synthase F0 subunit 6",
    "ATP8": "ATP synthase F0 subunit 8",
    "ND1": "NADH dehydrogenase subunit 1",
    "ND2": "NADH dehydrogenase subunit 2",
    "ND3": "NADH dehydrogenase subunit 3",
    "ND4": "NADH dehydrogenase subunit 4",
    "ND4L": "NADH dehydrogenase subunit 4L",
    "ND5": "NADH dehydrogenase subunit 5",
    "ND6": "NADH dehydrogenase subunit 6",
    "RNR1": "12S ribosomal RNA",
    "RNR2": "16S ribosomal RNA",
}

# EMMA-style tRNA /product strings ("tRNA-Trp(UCA)"), shared with rescue_trna.py
# via bin/mito_gene_order.py so both annotators emit byte-identical products.
# trna_product() falls back to a bare gene-name product if the suffix or anticodon
# isn't recognised -- a cosmetic field should never fail the run.
from mito_gene_order import TRNA_AA, trna_product


def map_gene_name(raw):
    """Return (emma_name, feature_type, frag, anticodon) for a MITOS2 BED
    feature name.

    feature_type is one of 'CDS', 'rRNA', 'tRNA'. ``frag`` is the MITOS
    fragment/exon part label ('a', 'b', '0', '1', ...) when the feature is one
    piece of a split gene, else None. ``anticodon`` is the tRNA anticodon in
    RNA form (e.g. 'UCA' from trnW(tca)), else None. Returns
    (None, None, None, None) for names we don't recognise (kept out of the
    EMMA contract rather than guessed).

    MITOS labels the exons of an intron-split gene with a trailing part suffix:
    a hyphen + letter (e.g. nad5-a, nad5-b) or an underscore + number
    (e.g. cox1_0, nad5_1). Anthozoan (coral) nad5 routinely carries a group I
    intron, so a split nad5 is normal, not an error. We strip the suffix so both
    parts map to the same EMMA gene and capture ``frag`` so the exons can later
    be concatenated in transcript order.
    """
    # Capture and strip the anticodon parens first (e.g. trnW(tca) -> trnW),
    # then peel off any MITOS fragment/exon suffix.
    anticodon_m = re.search(r"\(([a-zA-Z]+)\)$", raw.strip())
    anticodon = anticodon_m.group(1).upper().replace("T", "U") if anticodon_m else None
    bare_full = re.sub(r"\(.*?\)$", "", raw.strip())
    m = re.search(r"[-_](\d+|[a-z])$", bare_full)
    frag = m.group(1) if m else None
    bare = re.sub(r"[-_](?:\d+|[a-z])$", "", bare_full)

    if bare in PCG_MAP:
        return PCG_MAP[bare], "CDS", frag, None
    if bare in RRNA_MAP:
        return RRNA_MAP[bare], "rRNA", frag, None
    if bare.startswith("trn"):
        # trnF -> TF, trnL2 -> TL2, trnS1 -> TS1
        return "T" + bare[3:].upper(), "tRNA", frag, anticodon
    return None, None, None, None


def _frag_key(frag):
    """Sort key putting exons in MITOS transcript (part) order: a<b, 0<1."""
    if frag is None:
        return (0,)
    if frag.isdigit():
        return (int(frag),)
    return (ord(frag),)


def parse_bed(bed_path):
    """Parse MITOS2 result.bed, grouping the exons of split genes by gene name.

    Returns dict keyed by emma_name -> list of exon dicts
    {chrom, start, end, strand, ftype, frag} with 1-based inclusive start, the
    exon list ordered by MITOS transcript (part) order. A normal single-piece
    gene yields a one-element list; an intron-split gene (e.g. coral nad5) keeps
    one exon per BED line instead of being collapsed into a single span, so its
    reconstructed CDS skips the intron rather than swallowing it.
    """
    features = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                cols = line.split()
            if len(cols) < 6:
                continue
            chrom, start, end, raw_name, _score, strand = cols[:6]
            emma_name, ftype, frag, anticodon = map_gene_name(raw_name)
            if emma_name is None:
                continue
            start1, end = int(start) + 1, int(end)
            features.setdefault(emma_name, []).append({
                "chrom": chrom, "start": start1, "end": end,
                "strand": strand, "ftype": ftype, "frag": frag,
                "anticodon": anticodon,
            })
    for exons in features.values():
        exons.sort(key=lambda e: _frag_key(e["frag"]))
    return features


def _gff_span(start, end, seq_len):
    """Map a (possibly origin-spanning) feature to a GFF3-expressible span.

    GFF3 cannot express start > end on one line, so an origin-spanning feature
    (circular genome annotated by MITOS in circular mode) is reported up to the
    sequence end with the wrapped end recorded in a Note. Returns (gff_end, note).
    """
    if start > end:
        return seq_len, f";Note=origin-spanning;wrapped_end={end}"
    return end, ""


def write_gff(gff_path, features, chrom, seq_len, species, is_circular):
    # Order genes by their first (lowest) genomic coordinate.
    ordered = sorted(features.items(), key=lambda kv: min(e["start"] for e in kv[1]))
    with open(gff_path, "w") as out:
        out.write("##gff-version 3\n")
        out.write(f"##sequence-region\t{chrom}\t1\t{seq_len}\n")
        if species:
            out.write(f"##organism {species}\n")
        # Landmark region feature, matching EMMA's GFF. The Is_circular attribute
        # is what Geneious (and other GFF3 consumers) read to set the molecule
        # topology on import; without it a circular mitogenome comes in as linear.
        out.write(
            f"{chrom}\tmitos\tregion\t1\t{seq_len}\t.\t+\t.\t"
            f"ID={chrom}:1..{seq_len};Is_circular={'True' if is_circular else 'False'};"
            f"Name={chrom};genome=mitochondrion;mol_type=genomic DNA\n"
        )
        for emma_name, exons in ordered:
            product = PRODUCT.get(emma_name, emma_name)
            strand = exons[0]["strand"]
            ftype = exons[0]["ftype"]
            # Gene line spans the whole locus; the full (intron-skipping)
            # nucleotide/protein still goes to cds/ and proteins/ below, so
            # downstream BLAST/LCA and protein length stay correct.
            gene_start = min(e["start"] for e in exons)
            gene_end, note = _gff_span(gene_start, max(e["end"] for e in exons), seq_len)
            if len(exons) > 1:
                note = ";Note=multi-exon" + note
            attrs_gene = f"ID=gene-{emma_name};Name=MT-{emma_name};Product={product}{note}"
            out.write(f"{chrom}\tmitos\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attrs_gene}\n")
            # One feature line per exon, in transcript order. A single-exon gene
            # keeps the bare ID=<gene>; split genes get ID=<gene>.1, .2, ... so
            # each exon line is unique while sharing the same Parent.
            for i, e in enumerate(exons):
                feat_end, feat_note = _gff_span(e["start"], e["end"], seq_len)
                feat_id = emma_name if len(exons) == 1 else f"{emma_name}.{i + 1}"
                attrs_feat = (f"ID={feat_id};Parent=gene-{emma_name};"
                              f"Name=MT-{emma_name};Product={product}{feat_note}")
                out.write(f"{chrom}\tmitos\t{ftype}\t{e['start']}\t{feat_end}\t.\t{strand}\t.\t{attrs_feat}\n")


def _tbl_intervals(exons, seq_len):
    """Return the ordered list of (start, end) interval pairs for a .tbl
    feature block, in transcript order.

    Each exon's own start/end are plus-strand-referenced (as parse_bed and
    reorigin_features leave them; strand is tracked separately). start > end
    marks a feature that straddles the circular origin -- split into the two
    NCBI-tbl continuation intervals (start, seq_len) then (1, end), the same
    join syntax the format already uses for a multi-exon gene, so no separate
    Note is needed the way write_gff needs one for GFF3. A minus-strand
    interval is then written high..low (reversing the wrap-split order too),
    matching how EMMA's own .tbl encodes strand -- there is no strand column.

    ``exons`` must already be in MITOS transcript (part) order (parse_bed
    guarantees this), so multi-exon output here lists each exon's interval(s)
    in that same order without needing to re-sort.
    """
    pairs = []
    for e in exons:
        lo, hi, strand = e["start"], e["end"], e["strand"]
        plus_pts = [(lo, hi)] if lo <= hi else [(lo, seq_len), (1, hi)]
        pairs.extend([(b, a) for a, b in reversed(plus_pts)] if strand == "-" else plus_pts)
    return pairs


def write_tbl(tbl_path, features, chrom, seq_len, prefix, genetic_code):
    """Write an NCBI 5-column feature table (.tbl) for table2asn, from the same
    ``features`` dict write_gff consumes. Mirrors EMMA's .tbl shape (gene ->
    mRNA/CDS or gene -> tRNA/rRNA blocks) so process_files.py's normalisation
    (MT- stripping, protein_id removal, transl_except insertion) applies
    unchanged regardless of which annotator produced the table.
    """
    ordered = sorted(features.items(), key=lambda kv: min(e["start"] for e in kv[1]))
    with open(tbl_path, "w") as out:
        out.write(f">Feature {prefix}\n")
        for emma_name, exons in ordered:
            product = PRODUCT.get(emma_name, emma_name)
            ftype = exons[0]["ftype"]
            intervals = _tbl_intervals(exons, seq_len)

            def write_block(key, qual_lines):
                (s0, e0), rest = intervals[0], intervals[1:]
                out.write(f"{s0}\t{e0}\t{key}\n")
                for s, e in rest:
                    out.write(f"{s}\t{e}\n")
                for q in qual_lines:
                    out.write(f"\t\t\t{q[0]}\t{q[1]}\n")

            write_block("gene", [("gene", f"MT-{emma_name}")])
            if ftype == "CDS":
                write_block("mRNA", [("gene", f"MT-{emma_name}")])
                write_block("CDS", [("product", product), ("transl_table", genetic_code)])
            else:
                prod = trna_product(emma_name, exons[0].get("anticodon")) if ftype == "tRNA" else product
                write_block(ftype, [("product", prod)])


def gene_nt_seq(genome, chrom, start, end, strand):
    full = genome[chrom].seq
    # Origin-spanning feature (start > end): take the tail (start..end-of-seq)
    # followed by the head (1..end) so the gene is reconstructed across the
    # circular join before translation.
    seq = (full[start - 1:] + full[:end]) if start > end else full[start - 1:end]
    return seq.reverse_complement() if strand == "-" else seq


def anchor_origin(features, anchor):
    """Locate ``anchor``'s 5' end for re-origining, or None if MITOS did not call it.

    Returns (pos, rc) where ``pos`` is the 1-based genomic coordinate of the anchor's
    5' end in the current (pre-rotation) orientation and ``rc`` is True when the anchor
    sits on the minus strand, meaning the whole molecule must be reverse-complemented
    so the anchor reads 5'->3' on the plus strand from position 1 -- the orientation
    these records are deposited in.

    exons[0] is the 5'-most piece in MITOS transcript order, which matters for the
    multi-exon anchors: a tRNA is single-exon, but cox1 is intron-split in some sponge
    families and rrnL can be called in pieces too, and the origin belongs at the start
    of the gene, not the start of its largest exon.
    """
    exons = features.get(anchor)
    if not exons:
        return None
    e = exons[0]
    # On the minus strand the 5' end is the higher genomic coordinate (end).
    return (e["end"], True) if e["strand"] == "-" else (e["start"], False)


# Short, readable names for the log lines. The bare MITOS keys ("TM", "RNR2") are
# what the code and the asset speak, but an operator reading a run log should not
# have to translate them.
ANCHOR_LABEL = {
    "CO1": "cox1", "CO2": "cox2", "CO3": "cox3", "CYTB": "cytb",
    "RNR1": "rrnS (12S)", "RNR2": "rrnL (16S)",
}


def anchor_label(key):
    """Human-readable name for an anchor key, e.g. 'TM' -> 'tRNA-Met'.

    Not trna_product(): that needs an anticodon to produce 'tRNA-Met(CAU)' and
    degrades to the bare key ('tRNA-TM') without one, which is no more readable
    than the key itself. Here only the amino acid is wanted.
    """
    if key in ANCHOR_LABEL:
        return ANCHOR_LABEL[key]
    if key and key.startswith("T") and len(key) > 1:
        aa = TRNA_AA.get(key[1:])
        return f"tRNA-{aa}" if aa else key
    return key


def _flip(strand):
    return "-" if strand == "+" else "+"


def reorigin_seq(seq, pos, rc, seq_len):
    """Re-origin a circular sequence so ``pos`` (the anchor's 5' end) becomes base 1.

    Mirrors reorigin_features: on a minus-strand anchor the molecule is reverse-
    complemented first, then rotated so the anchor's 5' end is position 1. Gene-agnostic
    by construction -- it only ever sees a coordinate and a strand flag.
    """
    if rc:
        seq = seq.reverse_complement()
        pos = seq_len - pos + 1
    idx = (pos - 1) % seq_len
    return seq[idx:] + seq[:idx]


def reorigin_features(features, pos, rc, seq_len):
    """Shift every feature coordinate into the anchor-origined frame.

    ``pos`` is the anchor's 5' end in the current orientation. When ``rc`` each feature
    is first remapped onto the reverse complement (coords mirrored, strand
    flipped) so the whole transform stays consistent with reorigin_seq; then all
    coordinates are rotated so ``pos`` becomes position 1. A feature that ends up
    with start > end now straddles the new origin and is handled by the existing
    _gff_span / gene_nt_seq origin-spanning logic.
    """
    origin = seq_len - pos + 1 if rc else pos
    shifted = {}
    for name, exons in features.items():
        out = []
        for e in exons:
            s, en, st = e["start"], e["end"], e["strand"]
            if rc:
                s, en, st = seq_len - e["end"] + 1, seq_len - e["start"] + 1, _flip(st)
            ns = ((s - origin) % seq_len) + 1
            ne = ((en - origin) % seq_len) + 1
            out.append({**e, "start": ns, "end": ne, "strand": st})
        shifted[name] = out
    return shifted


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, type=Path, help="MITOS2 result.bed")
    ap.add_argument("--genome", required=True, type=Path, help="Assembly FASTA")
    ap.add_argument("--prefix", required=True, help="Output prefix (id.tech.date.assembler.mitos<ver>)")
    ap.add_argument("--outdir", required=True, type=Path, help="EMMA-style output dir")
    ap.add_argument("--code", type=int, required=True,
                    help="Mito genetic code (from meta.genetic_code); no default -- the "
                         "caller must supply it. It is written as the transl_table "
                         "qualifier, so a stale default would mislabel the annotation.")
    ap.add_argument("--origin-gene", required=True, metavar="KEY",
                    help="MITOS/EMMA feature key the published genome is re-origined to "
                         "(CO1, RNR2, TM, TF, ...), or NONE to skip re-origining. Resolved "
                         "per taxon by the caller from assets/taxonomy/mito_origin_anchors.json. "
                         "REQUIRED, like --code and for the same kind of reason: there are two "
                         "call sites (MITOS2 and CORAL_ANNOTATION_FIX, one per sample), and a "
                         "default here would let one of them be updated without the other, "
                         "silently publishing FIX and PASS corals of the same species on "
                         "different origins. A missing flag must fail loudly instead.")
    ap.add_argument("--origin-fallback", default="CO1", metavar="KEY",
                    help="Anchor to use when --origin-gene is not annotated, or NONE to leave "
                         "the genome on the pre-rotation origin. Defaults to CO1, which is "
                         "where rotate_to_cox1.py already put the molecule -- and re-origining "
                         "to MITOS's own cox1 call is more precise than the tblastn "
                         "back-extrapolation the pre-rotation used.")
    ap.add_argument("--species", default="", help="Species name for ##organism")
    ap.add_argument("--linear", action="store_true",
                    help="Mark the genome as linear (Is_circular=False) in the GFF "
                         "region feature. Default: circular, matching GetOrganelle/MitoHiFi.")
    args = ap.parse_args()

    if not args.bed.exists():
        sys.exit(f"BED not found: {args.bed}")
    if not args.genome.exists():
        sys.exit(f"Genome not found: {args.genome}")

    genome = SeqIO.to_dict(SeqIO.parse(str(args.genome), "fasta"))
    if not genome:
        sys.exit(f"No sequences read from {args.genome}")

    features = parse_bed(args.bed)
    if not features:
        print(f"WARNING: no recognised features parsed from {args.bed}", file=sys.stderr)

    # Anchor the GFF header on the contig MITOS annotated (fall back to first).
    chrom = next(iter(features.values()))[0]["chrom"] if features else next(iter(genome))
    if chrom not in genome:
        # MITOS may rename the contig; fall back to the single assembly record.
        chrom = next(iter(genome))
        for exons in features.values():
            for e in exons:
                e["chrom"] = chrom
    seq_len = len(genome[chrom].seq)

    # Re-origin the published genome to its taxon's deposition anchor, reusing MITOS's
    # own call for that gene. The upstream cox1 pre-rotation already lifted any
    # intron-split gene (the anthozoan nad5 group I intron) off the linearisation point
    # so MITOS annotated every gene in order; here we just shift the genome + all
    # coordinates so base 1 is the anchor's 5' end. Only for circular molecules --
    # rotating a linear contig would be wrong.
    #
    # The ladder never fails the run on a missing gene: anchor, then fallback, then
    # leave it where the pre-rotation put it. For the rrnL and trnF taxa the fallback
    # is a real improvement over the old behaviour rather than a degradation -- a sponge
    # whose RNR2 MITOS missed now lands on MITOS's own cox1 call, which is more precise
    # than the tblastn back-extrapolation that produced the pre-rotation origin.
    origin, used = None, None
    if not args.linear and (args.origin_gene or "").upper() != "NONE":
        for candidate in (args.origin_gene, args.origin_fallback):
            if not candidate or candidate.upper() == "NONE":
                continue
            origin = anchor_origin(features, candidate)
            if origin is not None:
                used = candidate
                break

    if origin is not None:
        pos, rc = origin
        if used != args.origin_gene:
            print(f"[mitos_to_emma] WARNING: no {anchor_label(args.origin_gene)} annotated; "
                  f"falling back to {anchor_label(used)}.", file=sys.stderr)
        features = reorigin_features(features, pos, rc, seq_len)
        rotated = reorigin_seq(genome[chrom].seq, pos, rc, seq_len)
        genome[chrom] = SeqRecord(rotated, id=genome[chrom].id,
                                  description=genome[chrom].description)
        print(f"[mitos_to_emma] re-origined {chrom} to {anchor_label(used)} 5' end "
              f"(was pos {pos}{', minus strand -> reverse-complemented' if rc else ''}); "
              f"position 1 is now {anchor_label(used)}.")
    elif not args.linear and (args.origin_gene or "").upper() != "NONE":
        print(f"[mitos_to_emma] WARNING: neither {anchor_label(args.origin_gene)} nor "
              f"{anchor_label(args.origin_fallback)} annotated; leaving genome on the "
              f"rotate_to_cox1 origin.", file=sys.stderr)

    args.outdir.mkdir(parents=True, exist_ok=True)
    cds_dir = args.outdir / "cds"
    prot_dir = args.outdir / "proteins"
    cds_dir.mkdir(exist_ok=True)
    prot_dir.mkdir(exist_ok=True)

    # Genome FASTA under the EMMA prefix
    SeqIO.write(list(genome.values()), str(args.outdir / f"{args.prefix}.fa"), "fasta")

    write_gff(args.outdir / f"{args.prefix}.gff", features, chrom, seq_len,
              args.species, is_circular=not args.linear)

    write_tbl(args.outdir / f"{args.prefix}.tbl", features, chrom, seq_len,
              args.prefix, args.code)

    for emma_name, exons in features.items():
        # Concatenate exons in MITOS transcript (part) order so an intron-split
        # gene yields its full spliced CDS with the intron removed; a single-exon
        # gene is unchanged.
        nt = Seq("")
        for e in exons:
            nt += gene_nt_seq(genome, e["chrom"], e["start"], e["end"], e["strand"])
        # EMMA contract: the gene identity is carried by the filename only; the
        # FASTA header is the bare prefix (no gene / MT- prefix), matching EMMA's
        # output so downstream modules that key off the header keep working.
        header = f">{args.prefix}"
        with open(cds_dir / f"{emma_name}.{args.prefix}.fa", "w") as out:
            out.write(f"{header}\n{nt}\n")
        if exons[0]["ftype"] == "CDS":
            aa = Seq(str(nt)).translate(table=args.code)
            with open(prot_dir / f"MT-{emma_name}.{args.prefix}.fa", "w") as out:
                out.write(f"{header}\n{aa}\n")

    print(f"[mitos_to_emma] wrote {len(features)} features to {args.outdir}")


if __name__ == "__main__":
    main()
