"""Unit tests for the taxonomy grading in bin/reference_divergence_check.py.

Tests classify_divergence directly (pure function, no GenBank / biopython needed),
focusing on the new CROSS_ORDER tier -- a reference from a different order, the
signal to switch to reference-free assembly -- and its graceful degradation to the
existing genus/family grading when a sample order is not supplied.
"""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "reference_divergence_check", ROOT / "bin" / "reference_divergence_check.py"
)
rdc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rdc)


class ClassifyDivergenceTests(unittest.TestCase):
    def tier(self, **kw):
        return rdc.classify_divergence(**kw)[0]

    def test_cross_order_detected_with_sample_order(self):
        # OG1422-like: Leptobrama (Carangiformes) vs Paraplagusia (Pleuronectiformes).
        self.assertEqual(
            self.tier(sample_species="Leptobrama muelleri", ref_org="Paraplagusia japonica",
                      ref_genus="Paraplagusia", ref_family="Cynoglossidae",
                      ref_order="Pleuronectiformes", sample_family="Leptobramidae",
                      sample_order="Carangiformes"),
            "CROSS_ORDER")

    def test_degrades_to_family_without_sample_order(self):
        self.assertEqual(
            self.tier(sample_species="Leptobrama muelleri", ref_org="Paraplagusia japonica",
                      ref_genus="Paraplagusia", ref_family="Cynoglossidae",
                      ref_order="Pleuronectiformes", sample_family="Leptobramidae"),
            "DIFFERENT_FAMILY")

    def test_degrades_to_non_congeneric_without_family_or_order(self):
        self.assertEqual(
            self.tier(sample_species="Oxyconger leptognathus", ref_org="Muraenesox cinereus",
                      ref_genus="Muraenesox", ref_family="Muraenesocidae", ref_order="Anguilliformes"),
            "NON_CONGENERIC")

    def test_same_order_is_confamilial_not_cross_order(self):
        self.assertEqual(
            self.tier(sample_species="Cephalopholis argus", ref_org="Epinephelus latifasciatus",
                      ref_genus="Epinephelus", ref_family="Serranidae", ref_order="Perciformes",
                      sample_family="Serranidae", sample_order="Perciformes"),
            "CONFAMILIAL")

    def test_congeneric_trusted(self):
        self.assertEqual(
            self.tier(sample_species="Plectropomus maculatus", ref_org="Plectropomus leopardus",
                      ref_genus="Plectropomus", ref_family="Serranidae", ref_order="Perciformes",
                      sample_order="Perciformes"),
            "CONGENERIC")

    def test_order_from_lineage(self):
        self.assertEqual(rdc.order_from_lineage(
            ["Eukaryota", "Chordata", "Actinopteri", "Pleuronectiformes", "Cynoglossidae"]),
            "Pleuronectiformes")
        self.assertEqual(rdc.order_from_lineage(["Eukaryota", "Chordata"]), "")


if __name__ == "__main__":
    unittest.main()
