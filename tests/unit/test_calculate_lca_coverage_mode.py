"""--coverage_mode: what calculateLCA.py's --min_coverage gate is measured against.

The script came from eDNA, where the query is a short ASV inside a longer reference,
so BLAST's query cover was a sensible gate. Here the query is a whole gene and most
references are barcode fragments, so an exact end-to-end match to a 600 bp COI record
covers under 40% of a 1.6 kb CO1 query. 'overlap' measures aligned length against the
shorter of query and subject instead.
"""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_calculate_lca_empty_output import load_calculate_lca


def blast_row(qlen, slen, length, qcovhsp, pident=99.5, name="Munidopsis recta"):
    """One 22-column outfmt-6 row in the order the BLAST module writes it."""
    cols = ["q1", "gi|1|gb|X1.1|", "1234", name, name, "Eukaryota", str(pident),
            str(length), str(qlen), str(slen), "1", "0", "0", "1", str(length), "1",
            str(length), f"{name} mitochondrion", "0.0", "500", str(qcovhsp), str(qcovhsp)]
    return "\t".join(cols) + "\n"


class CoverageForGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = load_calculate_lca()

    def test_qcovhsp_mode_returns_blast_value(self):
        self.assertEqual(self.m.coverage_for_gate("qcovhsp", 36.0, 293, 1019, 293), 36.0)

    def test_overlap_passes_fragment_matched_end_to_end(self):
        # 293 bp 12S fragment, fully aligned, against a 1019 bp query: 36% query cover.
        self.assertEqual(self.m.coverage_for_gate("overlap", 36.0, 293, 1019, 293), 100.0)

    def test_overlap_rejects_short_piece_of_long_reference(self):
        # 300 bp of a 16 kb mitogenome against a 1.6 kb query: overlap is 300/1600.
        self.assertAlmostEqual(self.m.coverage_for_gate("overlap", 19.0, 300, 1600, 16000), 18.75)

    def test_overlap_uses_query_when_query_is_shorter(self):
        # eDNA-shaped: 170 bp ASV fully aligned to a 16 kb reference.
        self.assertEqual(self.m.coverage_for_gate("overlap", 100.0, 170, 170, 16000), 100.0)


class ProcessBlastFileGateTests(unittest.TestCase):
    """The gate inside process_blast_file, without taxonomy databases."""

    @classmethod
    def setUpClass(cls):
        cls.m = load_calculate_lca()

    def _analyser(self):
        a = object.__new__(self.m.BLASTLCAAnalyzer)
        a.config = self.m.Config()
        a.logger = mock.MagicMock()
        a.corrector = mock.MagicMock()
        a.corrector.correct_line.side_effect = lambda line: line
        lineage = mock.MagicMock()
        lineage.species = "Munidopsis recta"
        a.assigner = mock.MagicMock()
        a.assigner.find_species_info.return_value = ("Munidopsis", "recta", "ncbi", lineage)
        a.assigner.lookup_taxid.return_value = None
        return a

    def _run(self, rows, mode, cover):
        with tempfile.TemporaryDirectory() as tmp:
            blast = Path(tmp) / "blast.tsv"
            blast.write_text("".join(rows))
            hits, _ = self._analyser().process_blast_file(
                blast, 98.0, cover, Path(tmp) / "missing.csv", coverage_mode=mode)
        return len(hits.get("q1", []))

    def test_default_mode_keeps_eDNA_behaviour(self):
        rows = [blast_row(qlen=1019, slen=293, length=293, qcovhsp=36)]
        self.assertEqual(self._run(rows, "qcovhsp", 40), 0)

    def test_overlap_mode_keeps_fragment_and_drops_short_piece(self):
        rows = [blast_row(qlen=1019, slen=293, length=293, qcovhsp=36),
                blast_row(qlen=1600, slen=16000, length=300, qcovhsp=19)]
        self.assertEqual(self._run(rows, "overlap", 80), 1)

    def test_reported_query_coverage_is_blast_value(self):
        rows = [blast_row(qlen=1019, slen=293, length=293, qcovhsp=36)]
        with tempfile.TemporaryDirectory() as tmp:
            blast = Path(tmp) / "blast.tsv"
            blast.write_text("".join(rows))
            hits, _ = self._analyser().process_blast_file(
                blast, 98.0, 80, Path(tmp) / "missing.csv", coverage_mode="overlap")
        # hit_info[12] is the query_coverage carried into the outputs
        self.assertEqual(hits["q1"][0][12], 36.0)


if __name__ == "__main__":
    unittest.main()
