#!/usr/bin/env python3
"""
Shared vertebrate mitochondrial gene-order and tRNA-naming tables.

Pure stdlib: no Biopython, no BLAST. Nextflow bind-mounts the whole bin/ dir onto
PATH, so a sibling import resolves via sys.path[0] (same pattern as
bin/orf_utils.py and bin/geo_loc_name_utils.py). This module must import cleanly
in the stdlib-only tylerpeirce/psycopg2:0.1 container.

REF_GENES was previously copy-pasted into annotation_stats.py,
emma_rescue_gate.py, trna_rescue_gate.py and rescue_trna.py, each carrying a
"keep this in sync" comment that had already drifted (three of them said "change
both"/"change all three" while there were four copies). The gates and the QC step
must agree byte-for-byte on what "present and in order" means, so the list lives
here and nowhere else.

Agreeing on the LIST turned out not to be enough: they also have to agree on how
a GFF is read into an ordered gene list, and they did not. parse_gff_attributes
had three copies and genes_by_coord had two byte-identical ones, while
annotation_stats.py kept a third, different implementation inline that deduped a
repeated gene by FIRST LINE SEEN rather than by lowest start. On a gene split
across the origin those two readings disagree, so the same GFF could be judged
in-order by a rescue gate and out-of-order by the QC step. Both readers now live
here, for the same reason REF_GENES does.
"""

# Reference gene order (tRNA, rRNA, CDS) -- the standard vertebrate set.
REF_GENES = [
    "TF", "RNR1", "TV", "RNR2", "TL2", "ND1", "TI", "TQ",
    "TM", "ND2", "TW", "TA", "TN", "TC", "TY", "CO1", "TS2",
    "TD", "CO2", "TK", "ATP8", "ATP6", "CO3", "TG", "ND3",
    "TR", "ND4L", "ND4", "TH", "TS1", "TL1", "ND5", "ND6",
    "TE", "CYTB", "TT", "TP",
]

TRNA_GENES = {g for g in REF_GENES if g.startswith("T")}   # the 22 mt tRNAs
RRNA_GENES = {"RNR1", "RNR2"}
PCG_GENES = {g for g in REF_GENES if g not in TRNA_GENES | RRNA_GENES}  # 13 PCGs

# EMMA-style tRNA gene suffix (the part after the leading 'T', e.g. "W", "S1",
# "L2") -> 3-letter amino acid, for tRNA /product strings ("tRNA-Trp(UCA)").
TRNA_AA = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L1": "Leu", "L2": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S1": "Ser", "S2": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}

# Canonical anticodon (DNA form, as tRNAscan-SE reports it) for each REF tRNA.
# The RNA form of these is what appears in the /product string, so this table and
# TRNA_PRODUCT below are two views of the same fact.
TRNA_ANTICODON = {
    "TF": "GAA", "TV": "TAC", "TL2": "TAA", "TL1": "TAG", "TI": "GAT",
    "TQ": "TTG", "TM": "CAT", "TW": "TCA", "TA": "TGC", "TN": "GTT",
    "TC": "GCA", "TY": "GTA", "TS2": "TGA", "TS1": "GCT", "TD": "GTC",
    "TK": "TTT", "TG": "TCC", "TR": "TCG", "TH": "GTG", "TE": "TTC",
    "TT": "TGT", "TP": "TGG",
}


def trna_product(emma_name, anticodon):
    """Return an EMMA-style tRNA /product string, e.g. 'tRNA-Trp(UCA)'.

    ``anticodon`` may be given in DNA (TCA) or RNA (UCA) form; it is transcribed
    to RNA for the product string. Falls back to a bare gene-name product if the
    suffix or anticodon isn't recognised -- a cosmetic field should never fail
    the run.
    """
    suffix = emma_name[1:] if emma_name.startswith("T") else emma_name
    aa = TRNA_AA.get(suffix)
    if aa and anticodon:
        return f"tRNA-{aa}({anticodon.upper().replace('T', 'U')})"
    return f"tRNA-{emma_name}"


# REF tRNA name -> canonical /product string, e.g. 'TW' -> 'tRNA-Trp(UCA)'.
TRNA_PRODUCT = {
    name: trna_product(name, ac) for name, ac in TRNA_ANTICODON.items()
}


# ----------------------------------------------------------------- GFF reading
#
# Shared by annotation_stats.py, emma_rescue_gate.py and trna_rescue_gate.py. See
# the module docstring: the gates and the QC step have to agree on how a GFF
# becomes an ordered gene list, not just on what the reference order is.


def parse_gff_attributes(attr_str):
    """Parse a GFF9 attribute column into a dict, ignoring malformed entries."""
    return dict(
        item.split("=", 1)
        for item in attr_str.strip().split(";")
        if "=" in item
    )


def gene_entries_by_coord(gff_path):
    """Return [(name, start, end, strand)] for the GFF's `gene` lines, by start.

    A gene written as more than one `gene` line -- which is how an origin-spanning
    feature is represented -- is kept ONCE, at its LOWEST start. Keeping the first
    line seen instead puts such a gene at whichever half the file happened to list
    first, which can make the order check fail on an annotation that is in fact
    correctly ordered; and a failed order check both fails the QC gate and silently
    suppresses the rescue gates.

    Lowest-start is deterministic where first-seen was not, so for the usual
    coordinate-sorted GFF nothing changes, and for a split gene the result no
    longer depends on file order. The returned start/end are the LOWER fragment's,
    which is what makes the ordering meaningful; a caller wanting the true span of
    an origin-spanning gene has to reconstruct it from the GFF itself.

    The `MT-` prefix EMMA writes is stripped, so names match REF_GENES.
    """
    best = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            attrs = parse_gff_attributes(parts[8])
            name = attrs.get("Name")
            if not name:
                continue
            gene = name.replace("MT-", "")
            start = int(parts[3])
            if gene not in best or start < best[gene][1]:
                best[gene] = (gene, start, int(parts[4]), parts[6])
    return sorted(best.values(), key=lambda entry: entry[1])


def genes_by_coord(gff_path):
    """Gene names present in the GFF, ordered by genomic start."""
    return [entry[0] for entry in gene_entries_by_coord(gff_path)]


# --------------------------------------------------- accepted order variants
#
# Accepted non-canonical gene orders, keyed by taxon. Each entry rewrites one
# CONTIGUOUS run of REF_GENES. Only add a rule with published evidence that the
# rearrangement is real -- an unexplained order must keep failing the gate.
#
# This is a CURATED TABLE, not a generic tolerance, and that is the design. A rule
# like "any single adjacent tRNA transposition is fine" would pass the exact case
# this gate exists to catch: an assembly whose apparent transposition is really a
# missed tRNA, leaving a tRNA-sized hole where the gene should have been called.
# Only a known taxon with a known variant passes; an unexpected rearrangement in
# an unexpected group still fails and gets human eyes on it.
#
# Entries are (rule_id, canonical_block, variant_block, evidence_accession).
#
# MAINTENANCE HAZARD, because most-specific-wins is NOT composition. Once a genus
# row exists for a genus inside a family that also has a row, the family rule stops
# applying to that genus entirely. Adding a genus row for some unrelated reason
# would silently remove the family's rule from every member of that genus. The
# rule for whoever edits this table: a new genus row inside a rule-carrying family
# must repeat that family's rules alongside its own.

ORDER_VARIANTS = {
    # Scarine parrotfishes: trnM and trnQ transposed relative to the canonical
    # IQM, giving IMQ. Confirmed against independent USNM-voucher GenBank records
    # that annotate the same coordinates.
    ("genus", "Hipposcarus"): [("scarine_imq", ("TQ", "TM"), ("TM", "TQ"), "PZ311005.1")],
    ("genus", "Chlorurus"):   [("scarine_imq", ("TQ", "TM"), ("TM", "TQ"), "PZ234023.1")],
    # trnD and trnS2 transposed between CO1 and CO2.
    ("genus", "Diploprion"):  [("diploprion_ds", ("TS2", "TD"), ("TD", "TS2"), "PZ244822.1")],
}

# Rules that are NOT in the table yet, kept here so the shape and the evidence
# question survive rather than being rediscovered.
#
# Both have strong INTERNAL evidence -- cross-assembler concordance, and a single
# block rewrite reconciling every affected assembly exactly -- but the design rule
# above asks for published evidence, and neither has an accession yet.
#
# ANGUILLIFORM_ND6_TE: ND6+trnE translocated from upstream of CYTB to between trnT
# and trnP. Keyed per FAMILY, never at ("order", "Anguilliformes"): anguilliforms
# in Synaphobranchidae and Nemichthyidae are byte-identical to REF_GENES and pass
# today, so this is not an order-level synapomorphy and an order key would grant
# licence across two families where the canonical order demonstrably holds. Before
# it ships, pull a CONGRID, NETTASTOMATID, COLOCONGRID or MURAENESOCID reference
# (a synaphobranchid or nemichthyid one confirms nothing) and ask it exactly one
# question: does it annotate trnE adjacent to the relocated ND6, or upstream of
# CYTB? If the latter, the rewrite is wrong in its trnE half and needs reshaping,
# not just an accession.
#
#   ANGUILLIFORM_ND6_TE = ("anguilliform_nd6_te",
#                          ("ND6", "TE", "CYTB", "TT"),
#                          ("CYTB", "TT", "ND6", "TE"),
#                          "TODO-accession")
#   ("family", "Congridae"):       [ANGUILLIFORM_ND6_TE],
#   ("family", "Nettastomatidae"): [ANGUILLIFORM_ND6_TE],
#   ("family", "Colocongridae"):   [ANGUILLIFORM_ND6_TE],
#   ("family", "Muraenesocidae"):  [ANGUILLIFORM_ND6_TE],
#   # Blachea is Colocongridae but resolves to blank family and order in the
#   # samplesheet, so the family rows cannot reach it. A narrow genus row with an
#   # EXPIRY CONDITION: delete it once the samplesheet taxonomy carries
#   # Colocongridae. Do not treat it as precedent for genus rows that exist only
#   # to dodge missing taxonomy, and note it is the live instance of the
#   # maintenance hazard above -- safe only because it names the same rule the
#   # family row would have given it.
#   ("genus", "Blachea"):          [ANGUILLIFORM_ND6_TE],
#
# macrourid_te: trnE ALONE translocated from between ND6 and CYTB to after trnP.
# NOT the anguilliform rule -- ND6 stays in place here. This is the weaker of the
# two despite affecting more samples: every affected assembly is one assembler and
# one read type, and every one carries a 61-68 bp unannotated hole exactly where
# canonical trnE sits. A systematic tRNA-calling miss would replicate across a
# family exactly as faithfully as biology would. If a Macrourinae reference
# annotates trnE in that gap, this row must NOT merge and those samples are an
# annotation defect instead.
#
#   ("family", "Macrouridae"): [
#       ("macrourid_te",
#        ("TE", "CYTB", "TT", "TP"),
#        ("CYTB", "TT", "TP", "TE"),
#        "TODO-accession")],

# Most specific rank wins OUTRIGHT. Ranks do not compose.
RANK_PRECEDENCE = ("genus", "family", "order", "class")

# Values that mean "the taxonomy did not resolve". All four occur in real
# samplesheets -- see isUnresolvedTaxon in subworkflows/local/prepare_samplesheet,
# which exists because unresolved family/order only warns while unresolved class
# aborts -- and none of them may ever match a rule key.
_UNRESOLVED = {"", "unknown", "na", "none", "dropped"}


def _normalise_rank_value(value):
    """Normalise a taxon value for key comparison, or '' if it is unresolved.

    Tolerates '', None and [] as well as the literal sentinels, because all of
    them reach here from real samplesheets.
    """
    if value is None or isinstance(value, (list, tuple, set, dict)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in _UNRESOLVED else text


def variant_rules_for(taxon):
    """Rules for the FIRST rank in RANK_PRECEDENCE that matches, else [].

    Most specific wins outright: a genus row shadows the family row for that
    genus entirely rather than adding to it. Comparison is case-insensitive on
    the stripped value.
    """
    taxon = taxon or {}
    lookup = {
        (rank, key.lower()): rules
        for (rank, key), rules in ORDER_VARIANTS.items()
    }
    for rank in RANK_PRECEDENCE:
        value = _normalise_rank_value(taxon.get(rank))
        if not value:
            continue
        rules = lookup.get((rank, value.lower()))
        if rules:
            return rules
    return []


def accepted_orders_for(taxon):
    """Every gene order this taxon may legitimately show, canonical FIRST.

    Returns [(order_list, rule_id_or_None), ...]. The canonical order is always
    present and always first, because a variant rule ADDS an accepted order, it
    does not replace one: the rearrangement is real in the clade but not universal
    within it, and a canonical member of a rule-carrying taxon must keep passing
    unchanged rather than being failed for NOT having the variant. Checking
    canonical first is also what keeps order_variant honest -- it stays "no" unless
    the variant is the thing that actually matched.

    Full order lists rather than a special-case comparison, so every consumer keeps
    its existing logic and simply compares against more than one list.

    Raises ValueError on a malformed table entry -- at import time via the
    self-check below, so a bad edit fails loudly instead of silently mis-gating a
    whole run.
    """
    canonical = (list(REF_GENES), None)
    rules = variant_rules_for(taxon)
    if not rules:
        return [canonical]
    variant, rule_id = _apply_rules(rules)
    if variant == list(REF_GENES):
        return [canonical]
    return [canonical, (variant, rule_id)]


def ref_order_for(taxon):
    """The single most specific accepted order: the variant if one applies, else
    canonical. For callers that need ONE ordering, such as picking the flanking
    search window for a rescued gene. Callers deciding whether an observed order
    is acceptable want accepted_orders_for instead.
    """
    orders = accepted_orders_for(taxon)
    return orders[-1]


def _apply_rules(rules):
    """Rewrite REF_GENES by every rule in turn; returns (order, joined_rule_ids)."""
    order = list(REF_GENES)
    applied = []
    for rule_id, canonical_block, variant_block in ((r[0], r[1], r[2]) for r in rules):
        canonical_block = list(canonical_block)
        variant_block = list(variant_block)
        if sorted(canonical_block) != sorted(variant_block):
            raise ValueError(
                f"ORDER_VARIANTS rule '{rule_id}': variant block {variant_block} is "
                f"not a permutation of canonical block {canonical_block}")
        start = _contiguous_index(order, canonical_block)
        if start is None:
            raise ValueError(
                f"ORDER_VARIANTS rule '{rule_id}': canonical block {canonical_block} "
                "is not a contiguous run of the reference order")
        order[start:start + len(canonical_block)] = variant_block
        applied.append(rule_id)
    return order, ";".join(applied)


def _contiguous_index(order, block):
    """Index at which `block` occurs as a contiguous run of `order`, else None."""
    n = len(block)
    for i in range(len(order) - n + 1):
        if order[i:i + n] == block:
            return i
    return None


def _validate_order_variants():
    """Fail at import on a malformed table rather than mis-gating a run."""
    for (rank, key), rules in ORDER_VARIANTS.items():
        if rank not in RANK_PRECEDENCE:
            raise ValueError(
                f"ORDER_VARIANTS key ('{rank}', '{key}'): rank is not one of "
                f"{RANK_PRECEDENCE}")
        if _normalise_rank_value(key) == "":
            raise ValueError(
                f"ORDER_VARIANTS key ('{rank}', '{key}'): an unresolved-taxon "
                "sentinel can never be a valid rule key")
        for rule in rules:
            if len(rule) != 4:
                raise ValueError(
                    f"ORDER_VARIANTS ('{rank}', '{key}'): expected "
                    "(rule_id, canonical_block, variant_block, accession)")
        ref_order_for({rank: key})


_validate_order_variants()


def matching_order_for(present, taxon):
    """The accepted order that `present` is in, or (None, None) if it is in none.

    `present` is a gene-name list in genomic order, as genes_by_coord returns.
    An order matches when `present` equals that order restricted to the genes
    actually present -- the same test annotation_stats.py applies, factored here so
    the QC step and the rescue gates cannot drift on it. Canonical is tried first,
    so the returned rule_id is None unless the variant is what actually matched.
    """
    present_set = set(present)
    for order, rule_id in accepted_orders_for(taxon):
        if present == [g for g in order if g in present_set]:
            return order, rule_id
    return None, None


def taxon_from_args(args):
    """Build the rank -> value dict from a parsed argparse namespace.

    Every consumer takes the same four options and builds the same dict; doing it
    once here is the point of this module. Unresolved values are normalised away
    inside variant_rules_for, so passing '' is safe.
    """
    return {
        "genus": getattr(args, "genus", "") or "",
        "family": getattr(args, "family", "") or "",
        "order": getattr(args, "taxon_order", "") or "",
        "class": getattr(args, "class_name", "") or "",
    }


def add_taxon_arguments(parser):
    """Add --genus/--family/--order/--class to an argparse parser.

    Taxonomy is consulted ONLY to look up a curated gene-order variant. It never
    affects which genes are expected, only which orderings of them are accepted.
    """
    parser.add_argument("--genus", default="",
                        help="genus, the first whitespace token of "
                             "meta.nominal_species_id. Used only to look up a "
                             "curated gene-order variant.")
    parser.add_argument("--family", default="", help="taxonomic family. See --genus.")
    parser.add_argument("--order", dest="taxon_order", default="",
                        help="taxonomic order. See --genus.")
    parser.add_argument("--class", dest="class_name", default="",
                        help="taxonomic class. See --genus.")
    return parser


# --------------------------------------------- LCA lineage (advisory cross-check)

# Ranks the cross-check considers, coarsest last so a caller can walk them.
LCA_RANKS = ("genus", "family", "order", "class")


def lca_lineage_from_combined(path):
    """Consensus lineage from an lca_combined.<prefix>.tsv, as {rank: value}.

    A rank gets a value only when EVERY non-empty row agrees on it; a rank the
    regions disagree about is reported as absent rather than resolved by majority,
    because the point of this file is to be a second OPINION and a split opinion is
    not one. 'dropped' and 'Unknown' -- what calculateLCA.py writes when it declines
    a rank or finds no lineage -- both count as absent.

    Stdlib only, and NEVER raises: a missing, unreadable or malformed file returns
    {}. This feeds an advisory column that holds nothing, so it must not be capable
    of failing a run.
    """
    absent = {"", "dropped", "unknown", "na", "n/a", "none", "null"}
    try:
        import csv
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except Exception:  # noqa: BLE001 - advisory only, never fail the run
        return {}

    lineage = {}
    for rank in LCA_RANKS:
        seen = set()
        for row in rows:
            value = (row.get(rank) or "").strip()
            if value and value.lower() not in absent:
                seen.add(value)
        if len(seen) == 1:
            lineage[rank] = seen.pop()
    return lineage


def order_variant_taxon_check(taxon, lineage):
    """Did a second taxonomy source agree with the rank a variant rule matched on?

    Returns 'agree', 'disagree:<lca_taxon>', 'unresolved', or 'no' when no rule
    fired. Purely a recorded flag: it holds nothing and changes no verdict.

    The variant table RELAXES the QC gate, so the obvious conservative design is
    to require two independent taxonomy sources to concur before a rule fires. That
    was rejected as a PRECONDITION, because the samples the table exists to unblock
    include ones whose samplesheet taxon is wrong or unresolved -- requiring
    agreement would ship the mechanism without releasing anything. Recording the
    disagreement instead makes a rule applied to a mislabelled sample queryable
    after the fact rather than invisible, which is the part that actually matters.
    """
    rules = variant_rules_for(taxon)
    if not rules:
        return "no"
    taxon = taxon or {}
    for rank in RANK_PRECEDENCE:
        value = _normalise_rank_value(taxon.get(rank))
        if not value:
            continue
        lookup = {(r, k.lower()) for (r, k) in ORDER_VARIANTS}
        if (rank, value.lower()) not in lookup:
            continue
        # This is the rank the rule matched on. Ask the LCA about that same rank.
        lca_value = (lineage or {}).get(rank, "")
        if not lca_value:
            return "unresolved"
        if lca_value.strip().lower() == value.lower():
            return "agree"
        return f"disagree:{lca_value.strip()}"
    return "unresolved"
