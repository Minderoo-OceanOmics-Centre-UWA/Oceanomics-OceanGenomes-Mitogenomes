"""Unit tests for bin/rescue_trna.py.

The tRNAscan-SE call is not exercised (that needs the BioContainer); the guard
logic in rescue_one() and the tabular parser are tested directly. Pure stdlib.
"""

import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "rescue_trna", ROOT / "bin" / "rescue_trna.py")
rt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rt)


def args_ns(**over):
    ns = types.SimpleNamespace(
        min_score=20.0, max_overlap=0, min_len=30, max_len=100,
        chrom="chr", ann_dir=None)
    for k, v in over.items():
        setattr(ns, k, v)
    return ns


# A minimal present-gene map: TS1's REF neighbours are TH (left) and TL1 (right).
GENES = {
    "TH":  (11000, 11070, "+"),
    "TL1": (11400, 11470, "+"),
    "ND5": (11480, 13000, "+"),
}
# Existing feature intervals (gene + child lines); nothing in the 11070..11400 gap.
FEATURES = [(11000, 11070), (11400, 11470), (11480, 13000)]


def hit(**over):
    h = dict(type="Ser", anticodon="GCT", lo=11150, hi=11220, strand="+",
             score=55.0, intron=0, note="")
    h.update(over)
    return h


class GuardTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        (self.tmp / "x.gff").write_text("##gff-version 3\n")
        (self.tmp / "x.tbl").write_text(">Feature chr\n")
        self.a = args_ns(ann_dir=self.tmp)

    def _run(self, hits, target="TS1"):
        return rt.rescue_one(target, hits, dict(GENES), list(FEATURES), self.a)

    def test_clean_hit_is_rescued_and_spliced(self):
        state, tgt, msg = self._run([hit()])
        self.assertEqual((state, tgt), ("RESCUED", "TS1"))
        gff = (self.tmp / "x.gff").read_text()
        tbl = (self.tmp / "x.tbl").read_text()
        self.assertIn("\ttRNA\t", gff)
        self.assertIn("Name=MT-TS1", gff)
        self.assertIn("product\ttRNA-Ser(GCU)", tbl)

    def test_wrong_anticodon_skips(self):
        # Ser isotype but the TS2 anticodon -> not TS1
        state, _, _ = self._run([hit(anticodon="TGA")])
        self.assertEqual(state, "SKIP")

    def test_low_score_skips(self):
        self.assertEqual(self._run([hit(score=12.0)])[0], "SKIP")

    def test_intron_skips(self):
        self.assertEqual(self._run([hit(intron=5)])[0], "SKIP")

    def test_hit_outside_gap_skips(self):
        # Inside the ND5 CDS, not the TH..TL1 gap
        self.assertEqual(self._run([hit(lo=12000, hi=12070)])[0], "SKIP")

    def test_overlap_with_existing_feature_skips(self):
        # Overlaps the TH gene by a base
        self.assertEqual(self._run([hit(lo=11040, hi=11110)])[0], "SKIP")

    def test_length_out_of_range_skips(self):
        self.assertEqual(self._run([hit(lo=11150, hi=11350)])[0], "SKIP")  # 201 bp

    def test_already_annotated_skips(self):
        genes = dict(GENES, TS1=(11150, 11220, "+"))
        state, _, msg = rt.rescue_one("TS1", [hit()], genes, list(FEATURES), self.a)
        self.assertEqual(state, "SKIP")
        self.assertIn("already annotated", msg)

    def test_ambiguous_two_comparable_hits_skips(self):
        h1 = hit(lo=11150, hi=11220, score=55.0)
        h2 = hit(lo=11250, hi=11320, score=54.0)   # different gap slot, ~equal score
        self.assertEqual(self._run([h1, h2])[0], "SKIP")

    def test_missing_neighbour_skips(self):
        genes = {"TH": (11000, 11070, "+")}   # no right neighbour present
        state, _, _ = rt.rescue_one("TS1", [hit()], genes, list(FEATURES), self.a)
        self.assertEqual(state, "SKIP")


class MultiTargetStateTests(unittest.TestCase):
    """Two targets in one invocation must not be placed on top of each other.

    rescue_one() used to append to the .gff/.tbl without recording the placement
    in the caller's `genes`/`features`, so the second target's --max-overlap
    guard could not see the tRNA the first had just written (and neighbour_gap
    still saw the pre-rescue gene set).
    """

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        (self.tmp / "x.gff").write_text("##gff-version 3\n")
        (self.tmp / "x.tbl").write_text(">Feature chr\n")
        self.a = args_ns(ann_dir=self.tmp)

    def test_second_target_cannot_land_on_the_first(self):
        # TW and TA are REF-order neighbours, both missing; the scan offers each
        # a hit at the SAME coordinates (the kind of duplicate a covariance model
        # produces for adjacent tRNAs).
        genes = {"ND2": (5000, 6000, "+"), "TN": (6200, 6270, "+"),
                 "CO1": (6300, 7800, "+")}
        features = [(5000, 6000), (6200, 6270), (6300, 7800)]
        span = dict(lo=6100, hi=6170, strand="+", score=55.0, intron=0, note="")
        hits = [dict(type="Trp", anticodon="TCA", **span),
                dict(type="Ala", anticodon="TGC", **span)]

        first = rt.rescue_one("TW", hits, genes, features, self.a)
        second = rt.rescue_one("TA", hits, genes, features, self.a)

        self.assertEqual(first[0], "RESCUED")
        self.assertEqual(second[0], "SKIP", second)
        self.assertIn("guards", second[2])
        # Exactly one tRNA feature reached the GFF.
        self.assertEqual((self.tmp / "x.gff").read_text().count("\ttRNA\t"), 1)

    def test_first_placement_narrows_the_next_gap(self):
        # After TW is placed, it becomes TA's left neighbour, so a TA hit that
        # sits upstream of TW (wrong side) is now out of gap.
        genes = {"ND2": (5000, 6000, "+"), "TN": (6400, 6470, "+")}
        features = [(5000, 6000), (6400, 6470)]
        rt.rescue_one("TW", [dict(type="Trp", anticodon="TCA", lo=6200, hi=6270,
                                  strand="+", score=55.0, intron=0, note="")],
                      genes, features, self.a)
        self.assertIn("TW", genes)
        state, _t, msg = rt.rescue_one(
            "TA", [dict(type="Ala", anticodon="TGC", lo=6050, hi=6120,
                        strand="+", score=55.0, intron=0, note="")],
            genes, features, self.a)
        self.assertEqual(state, "SKIP", msg)


class ParserTests(unittest.TestCase):
    TEXT = (
        "Sequence  tRNA Bounds  tRNA  Anti  Intron Bounds  Inf\n"
        "Name      tRNA#  Begin  End  Type  Codon  Begin  End  Score  Note\n"
        "--------  -----  -----  ---  ----  -----  -----  ---  -----  ----\n"
        "chr        1     11150  11220  Ser  GCT    0      0    55.3\n"
        "chr        2      9000   8930  Pro  TGG    0      0    41.0   pseudo\n"
        "chr        3      5000   5071  Leu  TAG   12     45    60.0\n"
    )

    def test_parse_skips_headers_reads_strand_note_intron(self):
        hits = rt.parse_trnascan(self.TEXT)
        self.assertEqual(len(hits), 3)
        self.assertEqual((hits[0]["type"], hits[0]["anticodon"]), ("Ser", "GCT"))
        self.assertEqual(hits[0]["strand"], "+")
        self.assertEqual(hits[1]["strand"], "-")          # 9000 -> 8930
        self.assertIn("pseudo", hits[1]["note"])
        self.assertEqual(hits[2]["intron"], 12)


if __name__ == "__main__":
    unittest.main()


class NeighbourGapVariantOrderTests(unittest.TestCase):
    """The flanking search window must follow the ACCEPTED gene order.

    A rescued tRNA is searched for between its neighbours. In a clade whose order
    is genuinely different those neighbours are different genes, so a
    canonical-only window brackets the wrong stretch of sequence and the rescue
    looks for the tRNA in the wrong place.
    """

    def setUp(self):
        import mito_gene_order as mgo
        self.mgo = mgo
        # Coordinates for every gene except TM, which is the one being "rescued".
        # Laid out in the SCARINE order (TI TM TQ ND2), so the variant window is
        # the genuinely correct one.
        self.variant_order, _ = mgo.ref_order_for({"genus": "Chlorurus"})
        self.genes = {}
        pos = 100
        for g in self.variant_order:
            if g != "TM":
                self.genes[g] = (pos, pos + 50)
            pos += 100

    def test_the_variant_order_gives_the_variant_neighbours(self):
        # In IMQ, TM sits between TI and TQ.
        gap = rt.neighbour_gap("TM", self.genes, self.variant_order)
        self.assertIsNotNone(gap)
        lo, hi = gap
        self.assertEqual(lo, max(self.genes["TI"]))
        self.assertEqual(hi, min(self.genes["TQ"]))

    def test_the_canonical_order_brackets_a_different_span(self):
        # In canonical IQM, TM sits between TQ and ND2 -- a different, and here
        # wrong, window. This is the defect, stated as a test.
        gap = rt.neighbour_gap("TM", self.genes, None)
        self.assertIsNotNone(gap)
        lo, hi = gap
        self.assertEqual(lo, max(self.genes["TQ"]))
        self.assertEqual(hi, min(self.genes["ND2"]))
        self.assertNotEqual(gap, rt.neighbour_gap("TM", self.genes, self.variant_order))

    def test_a_target_absent_from_the_order_yields_no_window(self):
        self.assertIsNone(rt.neighbour_gap("NOT_A_GENE", self.genes, None))
