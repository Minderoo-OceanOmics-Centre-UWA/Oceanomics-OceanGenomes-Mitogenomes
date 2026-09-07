#!/usr/bin/env python3
"""
Shared ORF / codon helpers for the annotation QC gate and the coral fixer.

Pure stdlib: no Biopython, no BLAST. Nextflow bind-mounts the whole bin/ dir
onto PATH, so a sibling import resolves via sys.path[0] (same pattern as
geo_loc_name_utils.py). This module must import cleanly in the gate's
tylerpeirce/psycopg2:0.1 container, which has no Biopython.

One source of truth for the per-code mitochondrial start/stop codon sets. The
partial tables in bin/process_files.py (`_MITO_STARTS_BY_CODE`,
`_MITO_STOPS_BY_CODE`) enumerate only code 2; this module enumerates every NCBI
mitochondrial translation table the pipeline can route (2, 4, 5, 9, 13, 14, plus
21/24/33), and process_files.py reads its `mito_start_codons` / `mito_stop_codons`
from here.

Sources: NCBI "The Genetic Codes" (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/
wprintgc.cgi), tables 2, 4, 5, 9, 13, 14, 21, 24, 33. Start sets are the codons
each table documents as initiators (the "M" / "i" positions); a CDS beginning on
one of these needs no transl_except and does not raise SEQ_FEAT.StartCodon.
"""

# --------------------------------------------------------------------------- #
# Per-code codon tables.
#
# Stops: every mitochondrial code reads TGA as Trp (never a stop). Only the
# vertebrate code (2) and the ascidian code (13) read AGA/AGG as stops; codes 4,
# 5, 9, 21, 24, 33 use TAA/TAG only (AGA/AGG there are Arg or Ser). Code 14 reads
# TAA as Tyr, so its only stop is TAG.
# --------------------------------------------------------------------------- #

STOP_CODONS_BY_CODE = {
    2: ("TAA", "TAG", "AGA", "AGG"),
    4: ("TAA", "TAG"),
    5: ("TAA", "TAG"),
    9: ("TAA", "TAG"),
    13: ("TAA", "TAG", "AGA", "AGG"),
    14: ("TAG",),
    21: ("TAA", "TAG"),
    24: ("TAA", "TAG"),
    33: ("TAA",),
}

START_CODONS_BY_CODE = {
    2: ("ATT", "ATC", "ATA", "ATG", "GTG"),
    4: ("TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"),
    5: ("TTG", "ATT", "ATC", "ATA", "ATG", "GTG"),
    9: ("ATG", "GTG"),
    13: ("TTG", "ATA", "ATG", "GTG"),
    14: ("ATG",),
    21: ("ATG", "GTG"),
    24: ("TTG", "CTG", "ATG", "GTG"),
    33: ("TTG", "CTG", "ATG", "GTG"),
}

SUPPORTED_CODES = frozenset(START_CODONS_BY_CODE) & frozenset(STOP_CODONS_BY_CODE)


# There is deliberately NO fallback for an unknown table. Every caller now resolves
# the code from meta.genetic_code and validates it, so an unrecognised value means a
# wiring bug, not an exotic organism -- and the previous silent default was wrong in
# both directions: it claimed the vertebrate code while actually returning the
# code-4/5/9 stop set (code 2 also stops on AGA/AGG). Failing loudly is the only way
# a wrong table cannot reach a submitted annotation.
def _require(code):
    try:
        code = int(code)
    except (TypeError, ValueError):
        raise ValueError(f"genetic code {code!r} is not an integer")
    if code not in SUPPORTED_CODES:
        raise ValueError(
            f"genetic code {code} is not a supported mitochondrial translation "
            f"table (supported: {sorted(SUPPORTED_CODES)})"
        )
    return code


def start_codons(code):
    return START_CODONS_BY_CODE[_require(code)]


def stop_codons(code):
    return STOP_CODONS_BY_CODE[_require(code)]


# --------------------------------------------------------------------------- #
# Translation (stdlib only).
# --------------------------------------------------------------------------- #

_BASE_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# Standard codon table (transl_table 1). The mitochondrial codes below differ
# from this only at a handful of codons, applied as per-code overrides.
_STANDARD = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Per-code deltas from the standard table (only codons that differ).
_CODE_OVERRIDES = {
    2: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
    4: {"TGA": "W"},
    5: {"AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W"},
    9: {"AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W"},
    13: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W", "GGA": "G"},
    14: {"AAA": "N", "AGA": "S", "AGG": "S", "TAA": "Y", "TGA": "W"},
    21: {"TGA": "W", "ATA": "M", "AAA": "N", "AGA": "S", "AGG": "S"},
    24: {"AGA": "S", "AGG": "K", "TGA": "W"},
    33: {"AGA": "S", "AGG": "K", "TAA": "Y", "TGA": "W"},
}


def _codon_table(code):
    tab = dict(_STANDARD)
    tab.update(_CODE_OVERRIDES.get(int(code), {}))
    return tab


def translate(nt, code):
    """Translate a nucleotide string in frame 0. Trailing 1-2 nt are ignored.
    An unknown codon (contains N or non-ACGT) yields 'X'."""
    tab = _codon_table(code)
    nt = nt.upper()
    out = []
    for i in range(0, len(nt) - 2, 3):
        out.append(tab.get(nt[i:i + 3], "X"))
    return "".join(out)


def revcomp(s):
    return "".join(_BASE_COMP.get(b, "N") for b in reversed(s.upper()))


# --------------------------------------------------------------------------- #
# CDS classification (used by annotation_qc_gate.py).
# --------------------------------------------------------------------------- #

def classify_cds(nt, code):
    """Inspect a spliced CDS nucleotide string under `code`.

    Returns dict:
        start_ok        first codon is a documented initiator for this table
        stop_ok         last full codon is a stop, OR the sequence length is not
                        a multiple of 3 (a poly-A-completed truncated stop, which
                        table2asn accepts via aa:TERM -- see process_files.py)
        internal_stops  number of stop codons before the last codon
        aa              the translation (terminal stop stripped)
        length_nt       len(nt)
    """
    nt = (nt or "").upper()
    starts = start_codons(code)
    stops = stop_codons(code)
    n = len(nt)

    start_ok = n >= 3 and nt[:3] in starts

    full_codons = n // 3
    last_codon = nt[(full_codons - 1) * 3: full_codons * 3] if full_codons else ""
    ends_on_stop = last_codon in stops
    # A CDS whose length is not a multiple of 3 is only an acceptable
    # poly-A-completed stop (table2asn's aa:TERM, see process_files.py) when the
    # trailing 1-2 nt are themselves the prefix of a stop codon for this table.
    # Any other remainder is a frameshift or a mis-called boundary -- exactly what
    # this gate exists to catch -- so it must NOT be waved through as truncated.
    truncated_stop = bool(n % 3) and any(
        s.startswith(nt[full_codons * 3:]) for s in stops)
    stop_ok = ends_on_stop or truncated_stop

    aa = translate(nt, code)
    body = aa[:-1] if aa.endswith("*") else aa
    internal_stops = body.count("*")

    return {
        "start_ok": bool(start_ok),
        "stop_ok": bool(stop_ok),
        "internal_stops": internal_stops,
        "aa": aa.rstrip("*"),
        "length_nt": n,
    }


# --------------------------------------------------------------------------- #
# ORF boundary refinement (used by coral_fix_bed.py; generalised from
# rescue_emma_pcg.py::refine_orf, which was hard-coded to transl_table 2).
# --------------------------------------------------------------------------- #

def _extend_orf(window, start, down_limit, stops, tab):
    """Translate in-frame from `start` (0-based) to the first stop, or to the
    last full codon before down_limit (poly-A convention). Returns
    (end, aa, poly_a) or None."""
    n = len(window)
    aa = []
    poly_a = False
    end = None
    for j in range(start, n - 2, 3):
        codon = window[j:j + 3]
        if codon in stops:
            end = j + 3
            break
        if j + 3 > down_limit:
            end = j
            poly_a = True
            break
        aa.append(tab.get(codon, "X"))
    if end is None or not aa:
        return None
    return end, "".join(aa), poly_a


def refine_orf(window, s_from, down_limit, code, up_codons=12):
    """Refine a boundary hint on a plus-oriented window to a clean CDS.

    window      : plus-oriented nucleotide string of the search window
    s_from      : 1-based position in `window` where the ORF is thought to start
    down_limit  : 1-based position past which the CDS must not extend (the
                  downstream gene's start, mapped into window coordinates)
    code        : NCBI mitochondrial translation table
    up_codons   : how many codons either side of s_from to search for an initiator

    The hint may be too early (a mis-called boundary upstream of the true ATG,
    as MITOS routinely does) or too late, so candidate initiators are searched
    both upstream and downstream of the hint, in-frame, within up_codons. The
    chosen start is the candidate whose translation to the first in-frame stop
    has no internal stop, preferring the one closest to the hint (tie-break:
    longer protein). If no candidate gives a clean ORF, the CDS is reported
    5'-partial from the window origin.

    Returns (cds_start, cds_end, aa, start_codon, poly_a) with cds_start/cds_end
    1-based inclusive in window coordinates, or None if no ORF could be built.
    `aa` never includes the terminal stop.
    """
    starts = start_codons(code)
    stops = stop_codons(code)
    tab = _codon_table(code)
    window = window.upper()
    frame0 = (s_from - 1) % 3
    hint = s_from - 1  # 0-based
    limit_nt = up_codons * 3

    candidates = []
    j = frame0
    while j + 3 <= len(window):
        if abs(j - hint) <= limit_nt and window[j:j + 3] in starts:
            ext = _extend_orf(window, j, down_limit, stops, tab)
            if ext is not None:
                end, aa, poly_a = ext
                codon = window[j:j + 3]
                # Rank: reaches a real encoded stop first; then canonical ATG
                # over an alternative initiator (GTG, then the rest); then the
                # initiator closest to the hint; then the longer protein. This
                # snaps a mis-called boundary to the true ATG even when a
                # permissive code (4/5) has an alternative start codon nearer.
                rank = 0 if codon == "ATG" else (1 if codon == "GTG" else 2)
                candidates.append((poly_a, rank, abs(j - hint), -len(aa), j, end, aa))
        j += 3

    if candidates:
        candidates.sort()
        poly_a, _rank, _dist, _neg_len, start, end, aa = candidates[0]
        return (start + 1, end, aa, window[start:start + 3], poly_a)

    # No clean initiator: 5'-partial from the first in-frame base.
    ext = _extend_orf(window, frame0, down_limit, stops, tab)
    if ext is None:
        return None
    end, aa, poly_a = ext
    return (frame0 + 1, end, aa, "partial", poly_a)
