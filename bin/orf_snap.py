#!/usr/bin/env python3
"""
Reference-free ORF boundary refinement for mitochondrial protein-coding genes.

This machinery used to live inside bin/coral_fix_bed.py, which meant the only way
to reach it was to be routed down the coral branch. It is not coral-specific: it
reads a MITOS BED and the assembly, and nudges each PCG's boundaries to the
nearest clean ORF under the sample's own genetic code. Nothing in it knows or
cares what lineage the sample is.

That mattered. In NOVA_260909_LA, 17 assemblies carried a genuinely broken CDS
(13 of them Porifera, the rest Malacostraca and Asteroidea) and not one was a
coral, so not one could reach the only repair code in the pipeline. Five
Demospongiae shared the identical defect -- ND4L starting on TTT and running to
a sense codon -- which is one systematic MITOS boundary miscall, not five bad
assemblies.

bin/coral_fix_bed.py now imports from here and keeps only what is genuinely
coral: BLAST transfer of rrnL, nad5, cox1 and mtMutS from an Anthozoa reference.
bin/invert_fix_bed.py wraps this module alone, for every invertebrate.
"""

import re

from orf_utils import classify_cds, refine_orf, revcomp



# EMMA-contract PCG bare names in the MITOS BED.
#
# nad5 IS in this set. It used to be excluded as a module constant because the
# coral reference-transfer pass rebuilds nad5 from a GenBank reference and a
# second opinion from the snap would contradict it. But that is a fact about the
# coral fixer, not about nad5, and for every non-coral invertebrate nothing else
# touches nad5 at all -- so excluding it here left the gene permanently
# unrepairable outside Anthozoa. The exclusion is now the caller's to declare,
# via snap_pcg_rows(skip=...), which coral_fix_bed.py does.
PCG_BARE = {
    "cox1", "cox2", "cox3", "cob", "atp6", "atp8",
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
}



# ---- BED helpers -----------------------------------------------------------

TRNA_RE = re.compile(r"^trn", re.I)

def parse_bed(path):
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cols = line.split("\t")
        if len(cols) < 6:
            cols = line.split()
        if len(cols) < 6:
            continue
        rows.append(cols[:6])  # chrom,start0,end,name,score,strand
    return rows

def bare_name(name):
    return re.sub(r"[-_](?:\d+|[a-z])$", "", re.sub(r"\(.*?\)$", "", name.strip()))


def cds_from_row(gseq, start0, end, strand):
    """Spliced CDS nucleotides for one BED row, handling the circular origin.

    MITOS writes an origin-spanning feature with start > end (e.g. `17684 19` on
    a 17.7 kb coral): the CDS is the tail of the sequence followed by the head.
    """
    n = len(gseq)
    seg = gseq[start0:end] if start0 < end else gseq[start0:] + gseq[:end]
    return revcomp(seg) if strand == "-" else seg.upper()


def _linear_span(start0, end, n):
    """A BED (start0, end) mapped into a doubled-genome coordinate space, where
    an origin-spanning feature becomes a plain interval (start0, end + n)."""
    return (start0, end) if start0 < end else (start0, end + n)


# A PCG that shares a canonical overlap with its neighbour (ATP8 into ATP6,
# ND4L into ND4) must be allowed to snap into that neighbour rather than being
# clipped at its first base, so the neighbour bound is relaxed by this many nt
# for those pairs. Anything more would let a boundary error walk into the
# neighbouring gene.

CANONICAL_OVERLAP = {"atp8": 12, "nad4l": 12}

# Per-gene absolute floor (nt) on how far a snap may move the CDS length. A flat
# floor of 60 nt is over a third of atp8 (~200 nt) and a fifth of nad3/nad4l, so
# short PCGs get a proportionate one instead.
LEN_TOL_FLOOR = {"atp8": 18, "nad4l": 18, "nad3": 24}
DEFAULT_LEN_TOL_FLOOR = 60



def snap_pcg_rows(rows, gseq, code, window=36, len_tol=0.15, skip=(),
                  lineage=None, circular=True):
    """Per-PCG ORF-boundary snap, assembly-only (no reference).

    For every single-exon PCG row whose reconstructed CDS is not a clean ORF
    under `code`, build a small window bounded by the neighbouring rows and call
    orf_utils.refine_orf. When it yields a clean ORF within `len_tol` of the
    original length, return the rewritten (start0, end) for that row.

    All coordinate work happens in a doubled-genome space so a feature spanning
    the circular origin (MITOS writes those with start > end) is an ordinary
    interval; results are mapped back, restoring the start > end form if the
    snapped CDS still wraps. Before this, every origin-spanning PCG was reported
    `no window` and silently went unchecked -- which in practice meant atp8 was
    never checked on any coral, because ROTATE_ORIGIN puts cox1 at position 1 and
    leaves atp8 straddling the join.

    `lineage` (taxonomic class or phylum) only widens the accepted initiators via
    orf_utils.LINEAGE_ALT_STARTS, so omitting it can make the snap attempt a
    repair it need not, never the reverse.

    `circular=False` (a linear contig) turns the doubled-genome space off: windows
    are clamped to the contig, so a gene cut by an assembly break is never snapped
    onto sequence from the other end.

    Bare names in `skip` are left alone: a gene rebuilt by a reference-transfer
    join owns its own coordinates, and snapping the stale single-exon row it
    replaces would only produce a contradictory second opinion.

    Returns (fixes, actions, unrepaired):
        fixes      = {row_index: (new_start0, new_end)}
        actions    = human-readable strings for the status file
        unrepaired = bare names of PCGs found broken that could NOT be repaired,
                     so the caller does not report the sample as FIXED
    """
    fixes = {}
    actions = []
    unrepaired = []

    n = len(gseq)
    gseq2 = (gseq + gseq).upper() if circular else gseq.upper()

    # Row spans in doubled coordinates, for neighbour clipping. Each row is also
    # considered one genome-length later, so a neighbour on the far side of the
    # origin still bounds a wrapped feature.
    spans = []
    for r in rows:
        try:
            lo, hi = _linear_span(int(r[1]), int(r[2]), n)
        except (ValueError, IndexError):
            continue
        spans.append((lo, hi))
        if circular:
            spans.append((lo + n, hi + n))

    by_bare = {}
    for i, r in enumerate(rows):
        by_bare.setdefault(bare_name(r[3]), []).append(i)

    for bare, idxs in by_bare.items():
        if bare not in PCG_BARE or len(idxs) != 1:
            continue  # only single-exon PCGs
        if bare in skip:
            continue  # owned by a reference-transfer pass (nad5-style join)
        i = idxs[0]
        r = rows[i]
        try:
            raw_s0, raw_e1 = int(r[1]), int(r[2])
        except (ValueError, IndexError):
            continue
        strand = r[5] if len(r) > 5 else "+"
        s0, e1 = _linear_span(raw_s0, raw_e1, n)

        nt = cds_from_row(gseq, raw_s0, raw_e1, strand)
        # Pass the lineage so a documented alternative initiator (decapod CO1 on
        # CGA, say) reads as a real start and the row is left alone. Without it
        # the snap would "repair" a boundary that was already correct.
        info = classify_cds(nt, code, lineage=lineage, gene=bare)
        if info["start_ok"] and info["stop_ok"] and not info["internal_stops"]:
            continue  # already a clean ORF

        def broken(reason):
            actions.append(f"PCG-skip {bare} ({reason})")
            unrepaired.append(bare)

        # Neighbour bounds in doubled coordinates, relaxed by the canonical
        # overlap this gene is allowed to share with its neighbour.
        slack = CANONICAL_OVERLAP.get(bare, 0)
        left_bound = max([sp[1] - slack for sp in spans if sp[1] <= s0] + [0])
        right_bound = min([sp[0] + slack for sp in spans if sp[0] >= e1] + [len(gseq2)])
        win_lo = max(left_bound, s0 - window)
        win_hi = min(right_bound, e1 + window)
        if win_lo >= win_hi:
            broken("no window")
            continue

        win = gseq2[win_lo:win_hi]
        if strand == "-":
            win = revcomp(win)
            s_from = (win_hi - e1) + 1
        else:
            s_from = (s0 - win_lo) + 1

        orf = refine_orf(win, s_from, len(win), code)
        if orf is None:
            broken("no clean ORF in window")
            continue
        cds_lo, cds_hi, aa, start_codon, _poly_a = orf
        old_len = e1 - s0
        new_len = cds_hi - cds_lo + 1
        if start_codon == "partial":
            broken(f"ORF not clean: start={start_codon}")
            continue
        floor = LEN_TOL_FLOOR.get(bare, DEFAULT_LEN_TOL_FLOOR)
        if abs(new_len - old_len) > max(floor, int(len_tol * old_len)):
            broken(f"len {old_len}->{new_len} outside tol")
            continue

        # Map window coords back to doubled-genome 0-based half-open, then onto
        # the real genome, restoring the start > end form for a wrapped CDS.
        if strand == "-":
            new_s0 = win_hi - cds_hi
            new_e1 = win_hi - cds_lo + 1
        else:
            new_s0 = win_lo + cds_lo - 1
            new_e1 = win_lo + cds_hi

        if new_s0 >= n:
            new_s0 -= n
            new_e1 -= n
        if new_e1 > n:
            new_e1 -= n          # wraps the origin: BED start > end, as MITOS writes it

        if (new_s0, new_e1) == (raw_s0, raw_e1):
            # The boundaries are already where refine_orf would put them, yet the
            # CDS did not classify clean -- an internal defect the snap can't fix.
            unrepaired.append(bare)
            actions.append(f"PCG-skip {bare} (no better boundary)")
            continue
        fixes[i] = (new_s0, new_e1)
        actions.append(
            f"PCG-fix {bare} {raw_s0}-{raw_e1}->{new_s0}-{new_e1} "
            f"({new_len}nt {len(aa)}aa)")

    return fixes, actions, unrepaired

