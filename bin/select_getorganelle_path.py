#!/usr/bin/env python3
"""Choose which GetOrganelle path FASTA becomes the assembly.

GetOrganelle writes one ``<prefix>.<label>.K<k>.{complete,scaffolds}.graph<G>.<P>
[.repeat_pattern<R>].path_sequence.fasta`` per path, where G numbers the disconnected
graph sets it found and P the alternative paths (isomers / repeat resolutions)
within one set. The modules used to copy ``*1.1.*.fasta``, i.e. whatever GetOrganelle
numbered graph 1. That numbering is not stable between runs. OG3000's reseed found
the same two graph sets twice (a 16,377 bp mitogenome and a 270 bp fragment): on
2026-09-16 the mitogenome was graph 1, on 2026-11-02 the fragment was, and the run
published a 155 bp "circular genome" from byte-identical inputs.

Rule: take path 1 of each graph set (so the isomer choice within a set is the one
GetOrganelle has always given us), then the set with the most sequence, then fewer
records, then the file name. With a single graph set this is exactly the old
``graph1.1`` file.

Every candidate is listed in the evidence TSV. No candidates (GetOrganelle produced
no path) writes nothing and exits 0, so the caller's empty placeholder still marks
"failed to assemble".

usage: select_getorganelle_path.py --dir OUTPUT_DIR --prefix PREFIX --name NAME
           --out assembly.fasta --evidence path_select.tsv
"""

import argparse
import re
import sys
from pathlib import Path

PATH_RE = re.compile(
    r"\.(?P<kind>complete|scaffolds)\.graph(?P<graph>\d+)\.(?P<path>\d+)"
    r"(?:\.repeat_pattern\d+)?\.path_sequence\.fasta$"
)


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
            cur.append(line)
    return ["".join(r) for r in recs]


def candidates(directory, prefix):
    """Path 1 of each graph set, as dicts. Other paths are ignored on purpose."""
    found = []
    for p in sorted(Path(directory).glob(f"{prefix}.*.path_sequence.fasta")):
        m = PATH_RE.search(p.name)
        if not m or int(m.group("path")) != 1:
            continue
        seqs = read_fasta(p)
        found.append({
            "file": p,
            "kind": m.group("kind"),
            "graph": int(m.group("graph")),
            "records": len(seqs),
            "bp": sum(len(s) for s in seqs),
        })
    return found


def choose(cands):
    usable = [c for c in cands if c["bp"] > 0]
    if not usable:
        return None
    return min(usable, key=lambda c: (-c["bp"], c["records"], c["file"].name))


def write_renamed(src, path, name):
    """Copy `src` line for line, renaming every header -- what `sed "/^>/s/.*/>NAME/"` did."""
    with open(src) as fin, open(path, "w") as out:
        for line in fin:
            out.write(f">{name}\n" if line.startswith(">") else line)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True, help="GetOrganelle output directory.")
    p.add_argument("--prefix", required=True, help="GetOrganelle --prefix, without the trailing dot.")
    p.add_argument("--name", required=True, help="FASTA header for every record written.")
    p.add_argument("--out", required=True)
    p.add_argument("--evidence", required=True)
    args = p.parse_args()

    cands = candidates(args.dir, args.prefix)
    best = choose(cands)

    with open(args.evidence, "w") as ev:
        ev.write("sample\tfile\tkind\tgraph\trecords\tbp\tchosen\n")
        for c in cands:
            ev.write(f"{args.name}\t{c['file'].name}\t{c['kind']}\t{c['graph']}\t{c['records']}\t"
                     f"{c['bp']}\t{'yes' if c is best else 'no'}\n")

    if best is None:
        print("No GetOrganelle path sequence found; nothing written.", file=sys.stderr)
        return 0
    write_renamed(best["file"], args.out, args.name)
    if len(cands) > 1:
        print(f"{len(cands)} GetOrganelle graph sets; chose graph {best['graph']} "
              f"({best['bp']} bp, {best['records']} record(s)).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
