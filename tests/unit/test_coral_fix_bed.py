"""Unit tests for the per-PCG ORF snap in bin/coral_fix_bed.py.

snap_pcg_rows() is pure Python (it only needs orf_utils), but the module imports
Biopython at the top for the reference-transfer half, so the import is stubbed
when Biopython is absent -- the same way the script runs in the MITOS2
BioContainer, where Biopython is present.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path load below.
sys.path.insert(0, str(ROOT / "bin"))

def _load_with_stubbed_biopython():
    """Load coral_fix_bed.py, standing in for Biopython if it is absent.

    The stubs are removed from sys.modules immediately afterwards: leaving them
    there would make `import Bio` succeed for every other test module in the same
    discovery run, so modules that legitimately skip without Biopython would
    instead fail on the first real Bio.* attribute they touch.
    """
    added = []
    try:
        import Bio  # noqa: F401
    except ImportError:
        bio = types.ModuleType("Bio")
        bio.SeqIO = None
        seq = types.ModuleType("Bio.Seq")
        seq.Seq = object
        sys.modules["Bio"], sys.modules["Bio.Seq"] = bio, seq
        added = ["Bio", "Bio.Seq"]
    try:
        spec = importlib.util.spec_from_file_location(
            "coral_fix_bed", ROOT / "bin" / "coral_fix_bed.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    finally:
        for name in added:
            sys.modules.pop(name, None)


cfb = _load_with_stubbed_biopython()

CODE = 4  # Coelenterate (Anthozoa)


def bed_row(name, start0, end, strand="+"):
    return ["chr", str(start0), str(end), name, "0.0", strand]


def clean_cds(n_codons, start="ATG", stop="TAA"):
    return start + "ACC" * n_codons + stop


class LinearSpanTests(unittest.TestCase):
    def test_plain_interval_unchanged(self):
        self.assertEqual(cfb._linear_span(100, 400, 1000), (100, 400))

    def test_origin_spanning_is_unrolled(self):
        # MITOS writes a wrapped feature as start > end, e.g. `17684 19`.
        self.assertEqual(cfb._linear_span(900, 20, 1000), (900, 1020))

    def test_cds_from_row_concatenates_across_the_origin(self):
        g = "AAAACCCCGGGGTTTT"          # 16 nt
        self.assertEqual(cfb.cds_from_row(g, 12, 4, "+"), "TTTTAAAA")


class SnapTests(unittest.TestCase):
    def test_clean_pcg_is_left_alone(self):
        body = clean_cds(40)
        genome = "GG" * 30 + body + "GG" * 30
        rows = [bed_row("cox1", 60, 60 + len(body))]
        fixes, actions, unrepaired = cfb.snap_pcg_rows(rows, genome, CODE)
        self.assertEqual(fixes, {})
        self.assertEqual(actions, [])
        self.assertEqual(unrepaired, [])

    def test_boundary_three_codons_early_is_snapped(self):
        # The MITOS nad1 failure mode: the called start sits a few codons before
        # the real ATG, so the CDS neither starts nor ends in frame.
        body = clean_cds(40)
        genome = "GG" * 30 + "AAAATACCC" + body + "GG" * 30
        true_start = 60 + 9
        rows = [bed_row("nad1", 60, true_start + len(body) - 9)]
        fixes, actions, unrepaired = cfb.snap_pcg_rows(rows, genome, CODE)
        self.assertEqual(unrepaired, [], actions)
        self.assertEqual(len(fixes), 1, actions)
        new_s0, new_e1 = fixes[0]
        self.assertEqual(genome[new_s0:new_e1], body)

    def test_origin_spanning_pcg_is_checked_and_snapped(self):
        """Regression: atp8 wrapping the origin used to report 'no window'.

        ROTATE_ORIGIN puts cox1 at position 1, which leaves atp8 straddling the
        join on most corals, so before the doubled-coordinate handling atp8 was
        never actually checked on any coral -- and OG2377 shipped an atp8 with no
        stop codon.
        """
        body = clean_cds(40)                       # 126 nt
        head, tail = body[100:], body[:100]
        genome = head + "GG" * 200 + tail          # body wraps the origin
        n = len(genome)
        # atp8 called 6 nt (2 codons) into its true start, and wrapped.
        rows = [bed_row("atp8", n - len(tail) + 6, len(head))]

        fixes, actions, unrepaired = cfb.snap_pcg_rows(rows, genome, CODE)
        self.assertNotIn("no window", " ".join(actions))
        self.assertEqual(unrepaired, [], actions)
        self.assertEqual(len(fixes), 1, actions)
        new_s0, new_e1 = fixes[0]
        self.assertGreater(new_s0, new_e1, "snapped CDS should still wrap")
        self.assertEqual(cfb.cds_from_row(genome, new_s0, new_e1, "+"), body)

    def test_unrepairable_pcg_is_reported_not_silently_skipped(self):
        # A truncated cox1 with no clean ORF anywhere in the window: the snap
        # cannot fix it, and the caller must hear about it so the sample stays
        # PARTIAL rather than being reported FIXED (OG2361).
        genome = "GG" * 40 + "ACC" * 60 + "GG" * 40
        rows = [bed_row("cox1", 80, 80 + 180)]
        fixes, actions, unrepaired = cfb.snap_pcg_rows(rows, genome, CODE)
        self.assertEqual(fixes, {})
        self.assertEqual(unrepaired, ["cox1"], actions)

    def test_short_pcg_length_floor_is_proportionate(self):
        # The flat 60 nt floor was a third of atp8; the per-gene floor is 18.
        self.assertLess(cfb.LEN_TOL_FLOOR["atp8"], cfb.DEFAULT_LEN_TOL_FLOOR)
        self.assertLess(cfb.LEN_TOL_FLOOR["nad4l"], cfb.DEFAULT_LEN_TOL_FLOOR)


if __name__ == "__main__":
    unittest.main()
