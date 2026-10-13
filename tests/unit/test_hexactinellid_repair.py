"""Glass-sponge MITOS repairs (bin/hexactinellid_repair.py via invert_fix_bed.py).

The real fixture is NOVA_260909_LA OG2650 (linear, 19.6 kb, code 4): MITOS dropped
nad2, started atp9 54 nt late and stopped atp8 short of its stop. The expected
BED is the hand repair, checked against Vazella/Oopsacas/Caulophacus RefSeq.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"
DATA = ROOT / "tests" / "test_data" / "hexactinellid"
sys.path.insert(0, str(BIN))

import hexactinellid_repair as hr  # noqa: E402
from orf_snap import parse_bed  # noqa: E402


def genome(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def pcgs(rows):
    return sorted((r[3], int(r[1]), int(r[2]), r[5]) for r in rows
                  if r[3].startswith(("nad", "cox", "cob", "atp")))


class OG2650Tests(unittest.TestCase):
    def test_repairs_match_the_reference_checked_hand_fix(self):
        rows = parse_bed(DATA / "OG2650.mitos.bed")
        g = genome(DATA / "OG2650.fasta")
        _actions, fixed, unrepaired = hr.repair(rows, g, 4, lineage="Hexactinellida")
        self.assertEqual(pcgs(rows), pcgs(parse_bed(DATA / "OG2650.expected.bed")))
        self.assertEqual({gene for gene, _ in fixed}, {"nad2", "atp8", "atp9"})
        self.assertEqual(unrepaired, [])

    def test_through_invert_fix_bed_with_lineage(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            proc = subprocess.run(
                [sys.executable, str(BIN / "invert_fix_bed.py"),
                 "--bed", str(DATA / "OG2650.mitos.bed"), "--genome", str(DATA / "OG2650.fasta"),
                 "--code", "4", "--lineage", "Hexactinellida", "--linear",
                 "--out-bed", str(tmp / "o.bed"), "--status", str(tmp / "s.tsv")],
                capture_output=True, text=True,
                env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(pcgs(parse_bed(tmp / "o.bed")),
                             pcgs(parse_bed(DATA / "OG2650.expected.bed")))
            self.assertTrue((tmp / "s.tsv").read_text().startswith("FIXED"))

    def test_other_lineages_are_untouched(self):
        """The stage is lineage-scoped: a demosponge with the same BED gets only
        the generic snap, never an added nad2."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            subprocess.run(
                [sys.executable, str(BIN / "invert_fix_bed.py"),
                 "--bed", str(DATA / "OG2650.mitos.bed"), "--genome", str(DATA / "OG2650.fasta"),
                 "--code", "4", "--lineage", "Demospongiae", "--linear",
                 "--out-bed", str(tmp / "o.bed"), "--status", str(tmp / "s.tsv")],
                capture_output=True, text=True,
                env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"}, check=True)
            self.assertNotIn("nad2", [r[3] for r in parse_bed(tmp / "o.bed")])


class RefusalTests(unittest.TestCase):
    """Each rule declines rather than guess."""

    def test_atp8_without_a_nearby_stop_is_left_unrepaired(self):
        g = "ATG" + "GCT" * 40 + "GCT" * 60          # no stop anywhere
        rows = [["c", "0", "60", "atp8", "0", "+"]]
        _a, fixed, unrep = hr.repair(rows, g, 4)
        self.assertEqual(fixed, [])
        self.assertIn("atp8", [u[0] for u in unrep])

    def test_atp9_start_far_from_reference_length_is_refused(self):
        # Upstream ATG exists but the resulting CDS would be ~2x reference length.
        g = "TAA" + "ATG" + "GCT" * 150 + "ACA" + "GCT" * 60 + "TAA"
        s = 3 + 3 + 450
        rows = [["c", str(s), str(len(g)), "atp9", "0", "+"]]
        _a, fixed, unrep = hr.repair(rows, g, 4)
        self.assertEqual(fixed, [])
        self.assertIn("atp9", [u[0] for u in unrep])

    def test_no_nad2_added_without_a_long_orf(self):
        cox3 = "ATG" + "GCT" * 100 + "TAA"
        spacer = ("GCTTAA" * 200)                     # stops every other codon
        nad5 = "ATG" + "GCT" * 300 + "TAA"
        g = cox3 + spacer + nad5
        rows = [["c", "0", str(len(cox3)), "cox3", "0", "+"],
                ["c", str(len(cox3) + len(spacer)), str(len(g)), "nad5", "0", "+"]]
        _a, fixed, unrep = hr.repair(rows, g, 4)
        self.assertNotIn("nad2", [r[3] for r in rows])
        self.assertIn("nad2", [u[0] for u in unrep])

    def test_nad2_search_does_not_take_nad5s_own_orf(self):
        """The first version picked nad5 itself (a 1.9 kb 'overlap')."""
        cox3 = "ATG" + "GCT" * 100 + "TAA"
        spacer = "GCTTAA" * 200
        nad5 = "ATG" + "GCT" * 600 + "TAA"
        g = cox3 + spacer + nad5
        rows = [["c", "0", str(len(cox3)), "cox3", "0", "+"],
                ["c", str(len(cox3) + len(spacer)), str(len(g)), "nad5", "0", "+"]]
        hr.repair(rows, g, 4)
        self.assertNotIn("nad2", [r[3] for r in rows])

    def test_partial_stop_only_when_the_frame_stays_open_to_the_next_feature(self):
        # nad1 runs to a T one nt before trnL2: end on the T (length % 3 == 1).
        nad1 = "ATG" + "GCT" * 100 + "GCT" + "T"
        trn = "GGGGGGGGGG"
        g = nad1 + "C" + trn   # C, not A: only the lone T can be a partial stop
        rows = [["c", "0", str(len(nad1) - 4), "nad1", "0", "+"],
                ["c", str(len(nad1) + 1), str(len(g)), "trnL2(taa)", "0", "+"]]
        _a, fixed, _u = hr.repair(rows, g, 4)
        self.assertEqual(int(rows[0][2]), len(nad1))
        self.assertIn("nad1", [f[0] for f in fixed])

    def test_lineage_gate(self):
        self.assertTrue(hr.is_hexactinellid("Hexactinellida"))
        self.assertTrue(hr.is_hexactinellid(" hexactinellida "))
        self.assertFalse(hr.is_hexactinellid("Demospongiae"))
        self.assertFalse(hr.is_hexactinellid(None))


if __name__ == "__main__":
    unittest.main()
