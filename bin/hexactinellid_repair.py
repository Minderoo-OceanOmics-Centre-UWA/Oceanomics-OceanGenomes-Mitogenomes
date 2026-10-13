"""Glass-sponge (Hexactinellida) MITOS repairs, assembly-only.

MITOS makes the same three calls on every glass sponge in NOVA_260909_LA (OG2629,
OG2636, OG2650, OG2937, OG2983), and the RefSeq glass sponges in assets/refdb
(Vazella pourtalesii NC_028054, Oopsacas minuta NC_027419, Caulophacus iocasicus
NC_072601, Aphrocallistes vastus NC_010769) contradict all three:

  nad2  dropped, or a short fragment kept, because it overlaps nad5. The overlap is
        real: Vazella 70 bp, Oopsacas 70 bp, Aphrocallistes 58 bp. The rebuilt
        genes align over the full reference length (78-85% aa identity).
  atp9  started 54 nt late on a non-start codon. Every reference atp9 is 240 bp
        (79 aa + stop) and opens on ATG; the assembly's ORF does too, from the
        first ATG after the upstream in-frame stop.
  atp8  stopped short of its stop codon; the frame runs on 15-87 nt to a TAA.

A fourth case is a CDS that runs into the next feature with no stop at all
(OG2983 nad1): mitochondrial transcripts complete a T or TA to TAA by
polyadenylation, so the CDS ends on the last in-frame T/TA before that feature.

Everything is read from the sample's own assembly. Nothing is copied from a
reference, so this can only ever annotate an ORF the sample carries. Each rule
declines (and says why) rather than guess.
"""

from orf_utils import classify_cds, revcomp, start_codons, stop_codons

ATP9_REF_NT = 240          # 79 aa + stop in all four references
ATP9_TOL = 0.10
ATP8_STOP_SEARCH_NT = 120
NAD2_MIN_CODONS = 300      # references 386-466 aa
NAD2_MIN_SPACER_LEAD = 300  # nad2 must open well before nad5 ...
NAD2_MAX_INTO_NAD5 = 400    # ... and may only run a short way into it
PARTIAL_STOP_GAP_NT = 6     # a partial stop sits right against the next feature


def is_hexactinellid(lineage):
    return bool(lineage) and lineage.strip().lower() == "hexactinellida"


def repair(rows, gseq, code, lineage=None):
    """Repair rows in place. Returns (actions, fixed, unrepaired), where fixed and
    unrepaired are lists of (gene, detail)."""
    n = len(gseq)
    g = gseq.upper()
    rcg = revcomp(g)
    stops = set(stop_codons(code))
    starts = set(start_codons(code))
    actions, fixed, unrepaired = [], [], []

    def seq_for(strand):
        return g if strand == "+" else rcg

    def to_gene(r):
        s, e = int(r[1]), int(r[2])
        return (s, e) if r[5] == "+" else (n - e, n - s)

    def to_bed(s, e, strand):
        return (s, e) if strand == "+" else (n - e, n - s)

    def find(name):
        hits = [r for r in rows if r[3] == name and int(r[1]) < int(r[2])]
        return hits[0] if len(hits) == 1 else None

    def clean(r, gene):
        nt = seq_for(r[5])[slice(*to_gene(r))]
        info = classify_cds(nt, code, lineage=lineage, gene=gene)
        return info["start_ok"] and info["stop_ok"] and not info["internal_stops"]

    def move(r, gene, s, e, why):
        old = f"{r[1]}-{r[2]}"
        b = to_bed(s, e, r[5])
        r[1], r[2] = str(b[0]), str(b[1])
        fixed.append((gene, f"{old} -> {b[0]}-{b[1]} ({why})"))
        actions.append(f"HEX-fix {gene} {old}->{b[0]}-{b[1]} ({why})")

    # atp8: extend to the first in-frame stop.
    r = find("atp8")
    if r and not clean(r, "atp8"):
        q = seq_for(r[5]); s, e = to_gene(r)
        if q[e - 3:e] not in stops:
            hit = next((p for p in range(e, min(e + ATP8_STOP_SEARCH_NT, n) - 2, 3)
                        if q[p:p + 3] in stops), None)
            if hit is None:
                unrepaired.append(("atp8", f"no in-frame stop within {ATP8_STOP_SEARCH_NT} nt"))
            else:
                move(r, "atp8", s, hit + 3, f"end extended {hit + 3 - e} nt to in-frame {q[hit:hit + 3]}")

    # atp9: first ATG after the upstream in-frame stop, at reference length.
    r = find("atp9")
    if r and not clean(r, "atp9"):
        q = seq_for(r[5]); s, e = to_gene(r)
        p = s
        while p - 3 >= 0 and q[p - 3:p] not in stops:
            p -= 3
        atg = next((x for x in range(p, e, 3) if q[x:x + 3] == "ATG"), None)
        if atg is None or atg == s:
            unrepaired.append(("atp9", "no upstream ATG in frame"))
        elif abs((e - atg) - ATP9_REF_NT) > ATP9_TOL * ATP9_REF_NT:
            unrepaired.append(("atp9", f"ATG start gives {e - atg} nt, reference is {ATP9_REF_NT}"))
        else:
            move(r, "atp9", atg, e, f"start moved {s - atg} nt upstream to ATG, {(e - atg) // 3 - 1} aa")

    # nad2: longest ORF opening in the cox3-nad5 spacer.
    c3, n5 = find("cox3"), find("nad5")
    old = find("nad2")
    if c3 and n5 and c3[5] == n5[5] and (old is None or not clean(old, "nad2")):
        strand = c3[5]; q = seq_for(strand)
        _, c3e = to_gene(c3); n5s, _ = to_gene(n5)
        best = None
        for fr in range(3):
            p = max(c3e - 30, 0) + fr
            while p < min(n5s + NAD2_MAX_INTO_NAD5, n) - 2:
                q0 = p
                while p < n - 2 and q[p:p + 3] not in stops:
                    p += 3
                end = p + 3
                if end <= n5s + NAD2_MAX_INTO_NAD5 and p < n - 2:
                    cands = [x for x in range(q0, p, 3)
                             if q[x:x + 3] in starts and c3e - 30 <= x < n5s - NAD2_MIN_SPACER_LEAD]
                    if cands:
                        # Prefer ATG, then another AT-start, close to the first
                        # candidate: code 4 also allows TTA/TTG/CTG/GTG, which would
                        # otherwise win just by sitting a few codons further upstream.
                        near = [x for x in cands if x - cands[0] < 30]
                        atg = [x for x in near if q[x:x + 3] == "ATG"]
                        atn = [x for x in near if q[x:x + 2] == "AT"]
                        s0 = (atg or atn or cands)[0]
                        if best is None or end - s0 > best[1] - best[0]:
                            best = (s0, end)
                p += 3
        if best and (best[1] - best[0]) // 3 - 1 >= NAD2_MIN_CODONS:
            if old is not None:
                rows.remove(old)
                actions.append(f"HEX-remove nad2 fragment {old[1]}-{old[2]}")
            b = to_bed(best[0], best[1], strand)
            rows.append([c3[0], str(b[0]), str(b[1]), "nad2", "0.0", strand])
            why = (f"ORF {q[best[0]:best[0] + 3]}..{q[best[1] - 3:best[1]]}, "
                   f"{(best[1] - best[0]) // 3 - 1} aa between cox3 and nad5, "
                   f"overlapping nad5 by {max(best[1] - n5s, 0)} nt")
            fixed.append(("nad2", f"added {b[0]}-{b[1]} ({why})"))
            actions.append(f"HEX-add nad2 {b[0]}-{b[1]} ({why})")
        else:
            unrepaired.append(("nad2", f"no ORF of >= {NAD2_MIN_CODONS} codons in the cox3-nad5 spacer"))

    # Partial stop: a CDS whose frame meets the next feature with no stop.
    feats = sorted((to_gene(x) + (x,)) for x in rows if int(x[1]) < int(x[2]))
    for r in [x for x in rows if x[3] in ("nad1", "nad2", "nad3", "nad4", "nad4l", "nad5",
                                         "nad6", "cob", "cox1", "cox2", "cox3", "atp6")]:
        if int(r[1]) >= int(r[2]) or clean(r, r[3]):
            continue
        q = seq_for(r[5]); s, e = to_gene(r)
        if q[e - 3:e] in stops:
            continue
        nxt = min((fs for fs, _fe, fr in feats if fr is not r and fr[5] == r[5] and fs >= e - 3),
                  default=None)
        if nxt is None:
            continue
        # The frame must stay open all the way to the next feature ...
        if any(q[p:p + 3] in stops for p in range(e, nxt - 2, 3)):
            continue
        # ... and end on T or TA within a few nt of it.
        for end in range(nxt, max(e, nxt - PARTIAL_STOP_GAP_NT) - 1, -1):
            rem = (end - s) % 3
            if (rem == 1 and q[end - 1] == "T") or (rem == 2 and q[end - 2:end] == "TA"):
                move(r, r[3], s, end, f"partial stop {q[end - rem:end]} {nxt - end} nt before the next feature")
                break

    rows.sort(key=lambda x: int(x[1]))
    return actions, fixed, unrepaired
