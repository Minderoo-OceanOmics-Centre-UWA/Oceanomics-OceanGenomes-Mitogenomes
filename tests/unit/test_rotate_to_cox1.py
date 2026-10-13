"""ROTATE_ORIGIN must never rotate a linear molecule.

A linear contig's two ends are assembly breaks, not neighbours. rotate_to_cox1.py
used to guard only on "more than one record", so every single linear contig
(including each SANITISE_FASTA _concat) was cut at cox1 and its ends glued into a
junction that does not exist. NOVA_260909_LA OG3043 reached MITOS2 rotated by
1300 bp even though MITOS2 was told --linear.
"""

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
SCRIPT = BIN / "rotate_to_cox1.py"

# Two arbitrary halves: if the script rotated at all, the output would start
# somewhere other than HALF_A.
HALF_A = "ATGACTAACTTATTCCCGGAAGGTACCATTAAAGCTGCC" * 5
HALF_B = "GGTTACCCAAGGTTTAACCCGGGTATATCGCGAATTCC" * 5


def run(seq, linear, n_records=1):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fa = "".join(f">c{i}\n{seq}\n" for i in range(n_records))
        (tmp / "g.fa").write_text(fa)
        (tmp / "cox1.faa").write_text(">cox1\nMTNLFPEGTIKAA\n")
        cmd = [sys.executable, str(SCRIPT), "--genome", str(tmp / "g.fa"),
               "--cox1-ref", str(tmp / "cox1.faa"), "--out", str(tmp / "out.fa")]
        if linear:
            cmd.append("--linear")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        out = (tmp / "out.fa").read_text() if (tmp / "out.fa").exists() else ""
        return proc, out


def seqs(fasta_text):
    out, cur = [], None
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            cur = []
            out.append(cur)
        elif cur is not None:
            cur.append(line.strip())
    return ["".join(c) for c in out]


class LinearPassthroughTests(unittest.TestCase):
    def test_linear_single_contig_is_written_through_unrotated(self):
        proc, out = run(HALF_A + HALF_B, linear=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(seqs(out), [HALF_A + HALF_B])
        self.assertIn("linear", proc.stdout.lower())

    def test_linear_flag_needs_no_blast(self):
        """The linear passthrough happens before any tblastn call, so it works
        even where BLAST is absent (and costs nothing where it is present)."""
        proc, _ = run(HALF_A + HALF_B, linear=True)
        self.assertNotIn("tblastn", (proc.stdout + proc.stderr).lower())

    def test_multi_record_passthrough_unchanged(self):
        proc, out = run(HALF_A, linear=False, n_records=2)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(seqs(out), [HALF_A, HALF_A])

    @unittest.skipUnless(shutil.which("tblastn"), "tblastn not on PATH")
    def test_circular_without_cox1_hit_still_runs(self):
        """Without --linear the cox1 search still runs; with no confident hit the
        record is written through, so the old circular path is intact."""
        proc, out = run(HALF_A + HALF_B, linear=False)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(seqs(out)), 1)


if __name__ == "__main__":
    unittest.main()
