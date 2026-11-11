"""bin/select_getorganelle_path.py picks the largest GetOrganelle graph set, not graph 1.

OG3000's reseed found the same two graph sets on 2026-09-16 and 2026-11-02 (a
16,262 bp mitogenome path and a 155 bp side graph) but numbered them the other way
round the second time, and the old `*1.1.*.fasta` copy published the 155 bp one.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "select_getorganelle_path.py"
PREFIX = "OG3000.ilmn.260909.getorg1770reseed"


def path_file(directory, graph, path, seqs, kind="complete", repeat=None):
    rp = f".repeat_pattern{repeat}" if repeat is not None else ""
    f = directory / f"{PREFIX}.animal_mt.K115.{kind}.graph{graph}.{path}{rp}.path_sequence.fasta"
    f.write_text("".join(f">scaffold_{i}\n{s}\n" for i, s in enumerate(seqs, 1)))
    return f


def run(directory):
    out, ev = directory / "out.fa", directory / "ev.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--dir", str(directory), "--prefix", PREFIX,
         "--name", PREFIX, "--out", str(out), "--evidence", str(ev)],
        capture_output=True, text=True)
    rows = [line.split("\t") for line in ev.read_text().splitlines()[1:]] if ev.exists() else []
    return proc, out, rows


class SelectGetorganellePathTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def chosen(self, rows):
        return [r[1] for r in rows if r[6] == "yes"]

    def test_largest_graph_wins_when_numbered_second(self):
        path_file(self.dir, 1, 1, ["A" * 155])
        big = path_file(self.dir, 2, 1, ["C" * 16262])
        proc, out, rows = run(self.dir)
        self.assertEqual(proc.returncode, 0)
        self.assertEqual(out.read_text(), f">{PREFIX}\n{'C' * 16262}\n")
        self.assertEqual(self.chosen(rows), [big.name])
        self.assertEqual(len(rows), 2)

    def test_largest_graph_wins_when_numbered_first(self):
        big = path_file(self.dir, 1, 1, ["C" * 16262])
        path_file(self.dir, 2, 1, ["A" * 155])
        _, out, rows = run(self.dir)
        self.assertEqual(self.chosen(rows), [big.name])
        self.assertIn("C" * 100, out.read_text())

    def test_single_graph_is_copied_byte_for_byte_with_renamed_header(self):
        # Wrapped lines must survive unchanged: this is what `cp` + `sed` produced.
        f = self.dir / f"{PREFIX}.animal_mt.K115.complete.graph1.1.path_sequence.fasta"
        f.write_text(">1+2-3+(circular)\nACGT\nACGT\nAC\n")
        _, out, _ = run(self.dir)
        self.assertEqual(out.read_text(), f">{PREFIX}\nACGT\nACGT\nAC\n")

    def test_only_path_one_of_each_graph_is_considered(self):
        # Isomers / repeat resolutions of one graph set: keep GetOrganelle's first.
        first = path_file(self.dir, 1, 1, ["G" * 100], repeat=1)
        path_file(self.dir, 1, 2, ["T" * 500], repeat=2)
        path_file(self.dir, 1, 11, ["T" * 900], repeat=3)
        _, out, rows = run(self.dir)
        self.assertEqual([r[1] for r in rows], [first.name])
        self.assertEqual(out.read_text(), f">{PREFIX}\n{'G' * 100}\n")

    def test_equal_size_prefers_fewer_records(self):
        path_file(self.dir, 1, 1, ["A" * 500, "C" * 500], kind="scaffolds")
        one = path_file(self.dir, 2, 1, ["G" * 1000], kind="scaffolds")
        _, _, rows = run(self.dir)
        self.assertEqual(self.chosen(rows), [one.name])

    def test_no_paths_writes_nothing_and_succeeds(self):
        proc, out, rows = run(self.dir)
        self.assertEqual(proc.returncode, 0)
        self.assertFalse(out.exists())
        self.assertEqual(rows, [])

    def test_other_prefixes_are_ignored(self):
        other = self.dir / "OG3000.ilmn.260909.getorg1770.animal_mt.K115.complete.graph1.1.path_sequence.fasta"
        other.write_text(">x\n" + "A" * 20000 + "\n")
        mine = path_file(self.dir, 1, 1, ["C" * 300])
        _, _, rows = run(self.dir)
        self.assertEqual(self.chosen(rows), [mine.name])


if __name__ == "__main__":
    unittest.main()
