"""Scaffolds are joined only across a direct, unique GetOrganelle graph edge.

The two shapes that matter come from NOVA_260909_LA:
  OG2629  scaffold 2 must follow scaffold 1 reverse-complemented; the graph has
          exactly one edge between those ends (plus a self-loop on the repeat
          unit). The old file-order concat glued it in forwards.
  OG2636  the scaffold end branches to two unselected nodes, one reaching each
          orientation of the partner. A unique wrong-orientation sequence
          overlap exists, so only the graph can say "don't join".
"""

import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
SCRIPT = BIN / "join_scaffolds_by_graph.py"
K = 7

rnd = random.Random(2629)


def rand(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def edge(i, seq, succ=()):
    """fastg records for both orientations of edge i. succ: {'+': [...], '-': [...]}
    with entries like (2, '-')."""
    return i, seq, dict(succ)


def write_graph(tmp, edges, scaffolds, gfa_segments):
    seqs = {i: s for i, s, _ in edges}
    lines = []
    for i, s, succ in edges:
        for o, body in (("+", s), ("-", rc(s))):
            head = f"EDGE_{i}_length_{len(s)}_cov_10.0{'' if o == '+' else chr(39)}"
            nxt = ",".join(
                f"EDGE_{j}_length_{len(seqs[j])}_cov_10.0{'' if jo == '+' else chr(39)}"
                for j, jo in succ.get(o, []))
            lines.append(f">{head}{':' + nxt if nxt else ''};\n{body}\n")
    fastg = tmp / "s.extended_K7.assembly_graph.fastg"
    fastg.write_text("".join(lines))
    gfa = tmp / "s.selected_graph.gfa"
    gfa.write_text("".join(f"S\t{name}\t{seq}\n" for name, seq in gfa_segments))
    fa = tmp / "in.fa"
    fa.write_text("".join(f">s\n{s}\n" for s in scaffolds))
    return fa, fastg, gfa


def run(tmp, fa, fastg, gfa):
    out, ev = tmp / "out.fa", tmp / "ev.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--fasta", str(fa), "--fastg", str(fastg),
         "--gfa", str(gfa), "--name", "s", "--out", str(out), "--evidence", str(ev)],
        capture_output=True, text=True)
    recs = [l for l in out.read_text().split(">") if l.strip()]
    seqs = ["".join(r.splitlines()[1:]) for r in recs]
    return proc, seqs, ev.read_text()


class GraphJoinTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        # A's right end overlaps B reverse-complemented by K bp (the OG2629 shape).
        self.A = rand(60)
        self.B = rc(self.A[-K:] + rand(50))

    def tearDown(self):
        self._tmp.cleanup()

    def test_unique_edge_joins_with_overlap_trimmed_and_reverse_complements(self):
        edges = [edge(1, self.A, {"+": [(2, "-")], "-": []}),
                 edge(2, self.B, {"+": [(1, "-"), (2, "+")], "-": [(2, "-")]})]  # B self-loop
        fa, fastg, gfa = write_graph(self.tmp, edges, [self.A, self.B],
                                     [("1", self.A), ("2", self.B)])
        proc, seqs, ev = run(self.tmp, fa, fastg, gfa)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(seqs, [self.A + rc(self.B)[K:]])
        self.assertIn("join", ev)
        self.assertNotIn("N", seqs[0])

    def test_branch_to_unselected_nodes_is_not_joined(self):
        """OG2636: the end reaches the partner in both orientations via other nodes."""
        C, D = rand(30), rand(30)
        edges = [edge(1, self.A, {"+": [(3, "+"), (4, "-")]}),
                 edge(2, self.B, {}),
                 edge(3, C, {"+": [(2, "+")]}),
                 edge(4, D, {"-": [(2, "-")]})]
        fa, fastg, gfa = write_graph(self.tmp, edges, [self.A, self.B],
                                     [("1", self.A), ("2", self.B)])
        proc, seqs, ev = run(self.tmp, fa, fastg, gfa)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(seqs, [self.A, self.B], "unproven scaffolds stay as written, in file order")
        self.assertNotIn("\njoin", ev)

    def test_edge_without_exact_overlap_is_refused(self):
        B_bad = rand(57)
        edges = [edge(1, self.A, {"+": [(2, "-")]}),
                 edge(2, B_bad, {"+": [(1, "-")]})]
        fa, fastg, gfa = write_graph(self.tmp, edges, [self.A, B_bad],
                                     [("1", self.A), ("2", B_bad)])
        _proc, seqs, ev = run(self.tmp, fa, fastg, gfa)
        self.assertEqual(seqs, [self.A, B_bad])
        self.assertIn("refused", ev)

    def test_non_reciprocal_link_is_not_joined(self):
        """A's end reaches only B, but B's end also reaches a third node."""
        C = rand(30)
        edges = [edge(1, self.A, {"+": [(2, "-")]}),
                 edge(2, self.B, {"+": [(1, "-"), (3, "+")]}),
                 edge(3, C, {})]
        fa, fastg, gfa = write_graph(self.tmp, edges, [self.A, self.B],
                                     [("1", self.A), ("2", self.B)])
        _proc, seqs, _ev = run(self.tmp, fa, fastg, gfa)
        self.assertEqual(seqs, [self.A, self.B])

    def test_single_scaffold_and_unreadable_graph_pass_through(self):
        fa = self.tmp / "one.fa"
        fa.write_text(f">s\n{self.A}\n")
        gfa = self.tmp / "g.gfa"
        gfa.write_text("")
        fastg = self.tmp / "s.extended_K7.assembly_graph.fastg"
        fastg.write_text(">not a fastg header\nACGT\n")
        _proc, seqs, ev = run(self.tmp, fa, fastg, gfa)
        self.assertEqual(seqs, [self.A])
        fa.write_text(f">s\n{self.A}\n>s\n{self.B}\n")
        proc, seqs, ev = run(self.tmp, fa, fastg, gfa)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(seqs, [self.A, self.B])
        self.assertIn("passthrough", ev)


if __name__ == "__main__":
    unittest.main()
