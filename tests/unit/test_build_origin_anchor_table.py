"""Unit tests for the reference-DB gene mapper behind the origin-anchor table.

This is the test that catches the bug that actually happened. A first cut of the
anchor tally read only the `gene` column, and tRNA rows in
assets/refdb/*/features.tsv routinely have an EMPTY gene column with only
"tRNA-Met" in `product`. Every tRNA silently vanished from the tally, and the
table reported Echinoidea as cox1 when the measured answer is tRNA-Phe (75.6%).
A second cut mapped tRNAs by first letter, which collapses Phe into Pro and
Thr/Trp/Tyr into each other -- same wrong answer, different route.

Both bugs are invisible downstream: the surviving rows still look plausible.
So the contract tested here is the strict one -- every FIRST feature in every
tracked database must map -- rather than a percentage.

Pure stdlib: no Biopython needed, matching bin/mito_gene_order.py itself.
"""

import collections
import csv
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from mito_gene_order import is_valid_anchor, refdb_gene_to_emma  # noqa: E402

REFDB = ROOT / "assets" / "refdb"
FEATURE_FILES = sorted(REFDB.glob("*/*_mito_refdb.features.tsv"))


class TrnaMappingTests(unittest.TestCase):
    """The amino-acid map must be a real map, not a first-letter shortcut."""

    def test_confusable_amino_acids_stay_distinct(self):
        # Phe/Pro both start P; Thr/Trp/Tyr all start T. Getting these wrong is
        # what decides the echinoderm answer.
        cases = {
            "tRNA-Phe": "TF", "tRNA-Pro": "TP",
            "tRNA-Thr": "TT", "tRNA-Trp": "TW", "tRNA-Tyr": "TY",
            "tRNA-Ala": "TA", "tRNA-Arg": "TR", "tRNA-Asn": "TN", "tRNA-Asp": "TD",
            "tRNA-Cys": "TC", "tRNA-Gln": "TQ", "tRNA-Glu": "TE", "tRNA-Gly": "TG",
            "tRNA-His": "TH", "tRNA-Ile": "TI", "tRNA-Lys": "TK", "tRNA-Met": "TM",
            "tRNA-Val": "TV",
        }
        for product, expected in cases.items():
            with self.subTest(product=product):
                # The empty gene column is the realistic case, and the one that broke.
                self.assertEqual(refdb_gene_to_emma("tRNA", "", product), expected)

    def test_tRNAs_map_from_an_empty_gene_column(self):
        self.assertEqual(refdb_gene_to_emma("tRNA", "", "tRNA-Met"), "TM")
        self.assertIsNotNone(refdb_gene_to_emma("tRNA", "", "tRNA-Phe"))

    def test_leu_and_ser_copy_numbers_are_kept_when_recorded(self):
        self.assertEqual(refdb_gene_to_emma("tRNA", "trnL2", "tRNA-Leu"), "TL2")
        self.assertEqual(refdb_gene_to_emma("tRNA", "trnS1", "tRNA-Ser"), "TS1")

    def test_bare_leu_and_ser_are_not_valid_anchors(self):
        # The DBs very often record only "tRNA-Leu". That is a real tally but it is
        # not addressable in a MITOS annotation, which always emits TL1/TL2.
        bare = refdb_gene_to_emma("tRNA", "", "tRNA-Leu")
        self.assertEqual(bare, "TL")
        self.assertFalse(is_valid_anchor(bare))
        self.assertFalse(is_valid_anchor(refdb_gene_to_emma("tRNA", "", "tRNA-Ser")))


class ProteinAndRrnaMappingTests(unittest.TestCase):
    def test_synonyms_across_eight_phyla_converge(self):
        for spelling in ("COX1", "cox1", "COI", "COXI", "MT-CO1",
                         "cytochrome c oxidase subunit I"):
            with self.subTest(spelling=spelling):
                self.assertEqual(refdb_gene_to_emma("CDS", spelling, ""), "CO1")
        for spelling in ("COX3", "COIII", "cytochrome c oxidase subunit III"):
            with self.subTest(spelling=spelling):
                self.assertEqual(refdb_gene_to_emma("CDS", spelling, ""), "CO3")
        for spelling in ("cob", "CYTB", "cytochrome b"):
            with self.subTest(spelling=spelling):
                self.assertEqual(refdb_gene_to_emma("CDS", spelling, ""), "CYTB")

    def test_large_and_small_rrna_synonyms(self):
        for spelling in ("rnl", "l-rRNA", "16S ribosomal RNA", "rrnL", "RRN16"):
            with self.subTest(spelling=spelling):
                self.assertEqual(refdb_gene_to_emma("rRNA", "", spelling), "RNR2")
        for spelling in ("rns", "s-rRNA", "12S ribosomal RNA", "rrnS"):
            with self.subTest(spelling=spelling):
                self.assertEqual(refdb_gene_to_emma("rRNA", "", spelling), "RNR1")

    def test_nadh_subunits_including_4L(self):
        self.assertEqual(refdb_gene_to_emma("CDS", "nad5", ""), "ND5")
        self.assertEqual(refdb_gene_to_emma("CDS", "ND4L", ""), "ND4L")
        self.assertEqual(refdb_gene_to_emma("CDS", "", "NADH dehydrogenase subunit 4L"), "ND4L")

    def test_genuinely_non_standard_genes_return_none(self):
        # Octocoral mutS, mussel F-ORF and homing endonucleases are real features in
        # these databases but are not part of the 37-gene set and must never become
        # an anchor. Returning None (and being counted) is the correct answer.
        self.assertIsNone(refdb_gene_to_emma("CDS", "mutS", "MutS-like protein"))
        self.assertIsNone(refdb_gene_to_emma("CDS", "F-ORF", "female open reading frame"))
        self.assertIsNone(refdb_gene_to_emma("CDS", "", "homing endonuclease"))


@unittest.skipUnless(FEATURE_FILES, "no reference databases checked out")
class TrackedDatabaseTests(unittest.TestCase):
    """Run the mapper over the real tracked databases, not synthetic input."""

    @staticmethod
    def _first_features():
        """(group, accession) -> the lowest-start feature's (type, gene, product)."""
        out = {}
        for path in FEATURE_FILES:
            group = path.parent.name
            lowest = {}
            with open(path) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    m = re.match(r"(\d+)-", (row.get("parts") or "").split(",")[0])
                    if not m:
                        continue
                    start = int(m.group(1))
                    acc = row["accession"]
                    if acc not in lowest or start < lowest[acc][0]:
                        lowest[acc] = (start, (row.get("type"), row.get("gene"),
                                               row.get("product")))
            for acc, (_start, triple) in lowest.items():
                out[(group, acc)] = triple
        return out

    def test_every_first_feature_maps(self):
        """The strict contract: an unmapped first feature silently skews the table."""
        firsts = self._first_features()
        self.assertGreater(len(firsts), 2000, "expected ~2100 tracked records")
        unmapped = collections.Counter()
        for (group, acc), triple in firsts.items():
            if refdb_gene_to_emma(*triple) is None:
                unmapped[(group,) + triple] += 1
        self.assertEqual(
            dict(unmapped), {},
            "these first features did not map, so their records dropped out of the "
            "anchor tally entirely")

    def test_scleractinia_still_tallies_as_trnM(self):
        """Guard the row whose loss would change ENA identity for submitted corals.

        Not a re-derivation of the whole table (test_origin_anchor_table.py does
        that); just a floor on the tally the mapper feeds, so a mapper regression
        that hid tRNAs again would fail here with a clear reason.
        """
        anthozoa = [t for (g, _a), t in self._first_features().items() if g == "anthozoa"]
        keys = collections.Counter(refdb_gene_to_emma(*t) for t in anthozoa)
        self.assertGreaterEqual(
            keys["TM"], 30,
            "tRNA-Met has nearly vanished from the anthozoan first-feature tally; the "
            "mapper is probably dropping tRNAs again")


if __name__ == "__main__":
    unittest.main()
