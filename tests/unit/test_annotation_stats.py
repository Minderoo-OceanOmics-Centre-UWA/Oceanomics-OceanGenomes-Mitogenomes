"""Unit tests for the reduced-tRNA classification in bin/annotation_stats.py.

Cnidaria (corals, anemones, jellyfish) and Porifera (sponges) both get the
relaxed completeness bar -- see InvertTaxonGroups in lib/ for why these two
groups are handled together while the rest of the invertebrate batch is not.
"""
import importlib.util
import sys
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# annotation_stats.py imports Bio.SeqIO for GFF/protein parsing elsewhere in the
# module, but the function under test here (has_reduced_trna_expectation) is a
# pure string classifier that never touches it. Biopython isn't installed in
# this local unit-test environment (it only runs in the module's BioContainer),
# so stub the import rather than pull in the real dependency just to satisfy it.
if "Bio" not in sys.modules:
    bio_stub = types.ModuleType("Bio")
    bio_seqio_stub = types.ModuleType("Bio.SeqIO")
    bio_stub.SeqIO = bio_seqio_stub
    sys.modules["Bio"] = bio_stub
    sys.modules["Bio.SeqIO"] = bio_seqio_stub

SPEC = importlib.util.spec_from_file_location(
    "annotation_stats", ROOT / "bin" / "annotation_stats.py"
)
astats = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(astats)


class ReducedTrnaExpectationTests(unittest.TestCase):
    def test_cnidaria_classes_are_reduced_trna(self):
        for class_name in ["Anthozoa", "Hydrozoa", "Scyphozoa", "Cnidaria"]:
            self.assertTrue(astats.has_reduced_trna_expectation(class_name))

    def test_porifera_classes_are_reduced_trna(self):
        for class_name in ["Demospongiae", "Calcarea", "Hexactinellida", "Porifera"]:
            self.assertTrue(astats.has_reduced_trna_expectation(class_name))

    def test_case_and_whitespace_insensitive(self):
        self.assertTrue(astats.has_reduced_trna_expectation("  porifera  "))
        self.assertTrue(astats.has_reduced_trna_expectation("ANTHOZOA"))

    def test_other_invert_phyla_are_not_reduced_trna(self):
        for class_name in ["Gastropoda", "Bivalvia", "Malacostraca", "Pycnogonida",
                            "Asteroidea", "Echinoidea", "Actinopteri", None, ""]:
            self.assertFalse(astats.has_reduced_trna_expectation(class_name))


if __name__ == "__main__":
    unittest.main()
