#!/usr/bin/env python3
"""Join GetOrganelle scaffolds only where the assembly graph proves the join.

GetOrganelle sometimes stops at >=2 scaffolds even though its own assembly graph
links them: typically one scaffold is a tandem-repeat unit (a self-loop in the
graph) whose copy number it cannot resolve, so it declines to write one path.
Downstream, SANITISE_FASTA then pastes the scaffolds together in file order. That
keeps every scaffold in the orientation GetOrganelle happened to write it, so a
scaffold that belongs reverse-complemented is glued in backwards. In
NOVA_260909_LA OG2629 that flipped a 2.7 kb scaffold, put ND6 on the wrong strand
and failed the annotation gate, although the graph has a direct edge
(1001709- -> 1000675-) giving the correct molecule.

Sequence overlap alone cannot stand in for the graph. OG2636 has an exact, unique
114 bp overlap in the WRONG orientation (a fold-back inverted repeat at its end),
while OG2629's true overlap is itself a repeat. Only the graph separates them: in
OG2636 the scaffold end branches to two unselected nodes, one reaching each
orientation of the partner, so no join is made.

Rule: scaffold end X joins scaffold end Y only when, in the full GetOrganelle
graph (extended_K*.assembly_graph.fastg),
  - X's successors, ignoring edges inside X's own scaffold (self-loops/repeats),
    are exactly {Y}, and
  - Y's successors back, ignoring edges inside Y's own scaffold, are exactly {X},
  - and the two ends share the exact k-mer overlap the graph implies.
The overlap is trimmed; nothing is ever inserted (no N spacers -- the
reference-guided join that used them was switched off for breaking genes).
Scaffolds without such evidence are left as they are, in file order, for
SANITISE_FASTA to concatenate exactly as before.

Fail-safe: any parsing problem writes the input through unchanged (exit 0).

usage: join_scaffolds_by_graph.py --fasta scaffolds.fa --fastg raw.fastg
           --gfa selected_graph.gfa [--gfa ...] --name NAME --out joined.fa
           --evidence graph_join.tsv
"""

import argparse
import re
import sys
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLIP = {"+": "-", "-": "+"}


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    recs, cur = [], None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            cur = []
            recs.append(cur)
        elif cur is not None:
            cur.append(line.upper())
    return ["".join(r) for r in recs]


def write_fasta(path, name, seqs):
    with open(path, "w") as out:
        for s in seqs:
            out.write(f">{name}\n")
            for i in range(0, len(s), 80):
                out.write(s[i:i + 80] + "\n")


EDGE = re.compile(r"EDGE_(\d+)_length_\d+_cov_[\d.]+('?)")


def parse_fastg(path):
    """Return (edge sequences in + orientation, successor map {(id, o): {(id, o)}})."""
    seqs, succ = {}, {}
    cur_id, cur_o, buf = None, None, []

    def flush():
        if cur_id is not None:
            s = "".join(buf).upper()
            seqs[cur_id] = s if cur_o == "+" else rc(s)

    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            flush()
            head = line[1:].rstrip(";")
            parts = head.split(":", 1)
            m = EDGE.match(parts[0])
            if not m:
                raise ValueError(f"unrecognised fastg header: {line[:80]}")
            cur_id, cur_o, buf = m.group(1), "-" if m.group(2) else "+", []
            nexts = set()
            if len(parts) > 1 and parts[1]:
                for tok in parts[1].split(","):
                    n = EDGE.match(tok)
                    if n:
                        nexts.add((n.group(1), "-" if n.group(2) else "+"))
            succ[(cur_id, cur_o)] = nexts
        else:
            buf.append(line)
    flush()
    return seqs, succ


def gfa_segments(paths):
    """{segment sequence: [edge ids in path order]} from selected_graph.gfa S lines."""
    segs = {}
    for p in paths:
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if f[0] == "S" and len(f) > 2:
                segs[f[2].upper()] = f[1].split("_")
    return segs


def orient_edge(scaf, edge_seq, at_start):
    """Orientation of `edge_seq` as it sits at the start (or end) of the scaffold."""
    probe = min(len(edge_seq), 200)
    if at_start:
        if scaf.startswith(edge_seq[:probe]):
            return "+"
        if scaf.startswith(rc(edge_seq)[:probe]):
            return "-"
    else:
        if scaf.endswith(edge_seq[-probe:]):
            return "+"
        if scaf.endswith(rc(edge_seq)[-probe:]):
            return "-"
    return None


def kmer_from(fastg_path, gfa_paths):
    m = re.search(r"extended_K(\d+)\.", Path(fastg_path).name)
    if m:
        return int(m.group(1))
    for p in gfa_paths:
        for line in open(p):
            f = line.split("\t")
            if f[0] == "L" and len(f) > 5:
                mm = re.match(r"(\d+)M", f[5])
                if mm:
                    return int(mm.group(1))
    return None


def plan_joins(scafs, seqs, succ, segs, k):
    """Return (links, ends, notes). links: {(i, side): (j, side, k)} both directions."""
    notes = []
    ends = {}            # (i, 'L'|'R') -> oriented edge pointing OUT of the scaffold
    owned = {}           # i -> set of edge ids inside scaffold i
    for i, s in enumerate(scafs):
        path = segs.get(s) or segs.get(rc(s))
        if path is None:
            notes.append(f"scaffold {i + 1}: no matching selected_graph segment")
            continue
        first, last = path[0], path[-1]
        if first not in seqs or last not in seqs:
            notes.append(f"scaffold {i + 1}: edge missing from fastg")
            continue
        of = orient_edge(s, seqs[first], at_start=True)
        ol = orient_edge(s, seqs[last], at_start=False)
        if of is None or ol is None:
            # The segment may be written reverse-complemented relative to the path.
            first, last = path[-1], path[0]
            of = orient_edge(s, seqs[first], at_start=True)
            ol = orient_edge(s, seqs[last], at_start=False)
        if of is None or ol is None:
            notes.append(f"scaffold {i + 1}: could not orient its end edges")
            continue
        ends[(i, "R")] = (last, ol)               # leaving the right end
        ends[(i, "L")] = (first, FLIP[of])        # leaving the left end (reverse walk)
        owned[i] = set(path)

    # Which scaffold end does an oriented edge ENTER? Entering i's left end means
    # walking onto its first edge in its own orientation; entering its right end
    # means walking onto its last edge reversed.
    enters = {}
    for (i, side), (e, o) in ends.items():
        enters[(e, FLIP[o])] = (i, side)

    def partners(i, side):
        e, o = ends[(i, side)]
        out = {n for n in succ.get((e, o), set()) if n[0] not in owned[i]}
        return out

    links = {}
    for (i, side) in list(ends):
        out = partners(i, side)
        if len(out) != 1:
            continue
        hit = enters.get(next(iter(out)))
        if hit is None or hit[0] == i:
            continue
        j, jside = hit
        if len(partners(j, jside)) != 1:
            continue
        back = enters.get(next(iter(partners(j, jside))))
        if back != (i, side):
            continue
        links[(i, side)] = (j, jside)
    return links, ends, notes


def oriented(s, flip):
    return rc(s) if flip else s


def build(scafs, links, k):
    """Walk the unique links into chains; return (sequences, evidence rows)."""
    used, out, rows = set(), [], []
    order = list(range(len(scafs)))

    def walk_from(start):
        # Go left first to find the chain's true start (guard against cycles).
        i, flip, seen = start, False, {start}
        while True:
            side = "R" if flip else "L"
            nxt = links.get((i, side))
            if nxt is None or nxt[0] in seen:
                break
            j, jside = nxt
            # Arriving at j through jside, travelling leftwards: j's jside must be
            # its right end in the chain orientation.
            flip = (jside == "L")
            i = j
            seen.add(j)
        return i, flip

    for start in order:
        if start in used:
            continue
        i, flip = walk_from(start)
        seq = oriented(scafs[i], flip)
        chain = [f"{i + 1}{'-' if flip else '+'}"]
        used.add(i)
        while True:
            side = "L" if flip else "R"
            nxt = links.get((i, side))
            if nxt is None or nxt[0] in used:
                if nxt is not None:
                    rows.append(("cycle", chain[-1], f"{nxt[0] + 1}", "graph links close a circle; left open here"))
                break
            j, jside = nxt
            jflip = (jside == "R")
            piece = oriented(scafs[j], jflip)
            if seq[-k:] != piece[:k]:
                rows.append(("refused", chain[-1], f"{j + 1}{'-' if jflip else '+'}",
                             f"graph edge present but the {k} bp overlap is not exact"))
                break
            seq += piece[k:]
            chain.append(f"{j + 1}{'-' if jflip else '+'}")
            rows.append(("join", chain[-2], chain[-1], f"direct unique graph edge, {k} bp overlap trimmed"))
            used.add(j)
            i, flip = j, jflip
        out.append((seq, chain))
    return out, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--fastg", required=True, type=Path, help="raw extended_K*.assembly_graph.fastg")
    ap.add_argument("--gfa", required=True, type=Path, action="append",
                    help="selected_graph.gfa (repeatable)")
    ap.add_argument("--name", required=True, help="FASTA header for every output record")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--evidence", required=True, type=Path)
    args = ap.parse_args()

    scafs = read_fasta(args.fasta)
    header = "action\tleft\tright\tnote\n"

    def passthrough(reason):
        write_fasta(args.out, args.name, scafs)
        args.evidence.write_text(header + f"passthrough\t-\t-\t{reason}\n")
        print(f"[join_scaffolds_by_graph] passthrough: {reason}")
        return 0

    if len(scafs) < 2:
        return passthrough(f"{len(scafs)} scaffold(s): nothing to join")
    try:
        seqs, succ = parse_fastg(args.fastg)
        segs = gfa_segments(args.gfa)
        k = kmer_from(args.fastg, args.gfa)
        if not k:
            return passthrough("could not determine the graph k-mer")
        links, _ends, notes = plan_joins(scafs, seqs, succ, segs, k)
    except (OSError, ValueError) as exc:
        return passthrough(f"graph unreadable: {exc}")

    chains, rows = build(scafs, links, k)
    write_fasta(args.out, args.name, [s for s, _ in chains])
    with open(args.evidence, "w") as out:
        out.write(header)
        for r in rows:
            out.write("\t".join(r) + "\n")
        for n in notes:
            out.write(f"note\t-\t-\t{n}\n")
        for s, chain in chains:
            out.write(f"molecule\t{'>'.join(chain)}\t-\t{len(s)} bp\n")
    print(f"[join_scaffolds_by_graph] {len(scafs)} scaffold(s) -> {len(chains)} "
          f"({sum(r[0] == 'join' for r in rows)} graph join(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
