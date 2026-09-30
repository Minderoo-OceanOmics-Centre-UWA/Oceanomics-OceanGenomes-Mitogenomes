"""An _rgj / _concat is a scaffolding artefact and must never be submittable.

Nothing in the pipeline distinguished these until now: the QC gate is content-based and
submission_ready is pure gate arithmetic, so OG2936's _rgj reached submission_ready=true
in NOVA_260909_LA and would have been submitted to PRJEB110568.
"""

import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

from assembly_variant import RESCUE_SUFFIXES, is_rescue_molecule


class RescueMoleculeTests(unittest.TestCase):
    def test_bare_codes(self):
        self.assertTrue(is_rescue_molecule("getorg1770_rgj"))
        self.assertTrue(is_rescue_molecule("getorg1770_concat"))
        self.assertTrue(is_rescue_molecule("getorg1770reseed_rgj"))
        self.assertTrue(is_rescue_molecule("getorg1770reseed_concat"))

    def test_full_prefixes_and_seqids(self):
        """Callers hold different amounts of the name, so all three forms must work."""
        self.assertTrue(is_rescue_molecule("OG2936.ilmn.260909.getorg1770_rgj"))
        self.assertTrue(is_rescue_molecule("OG2936.ilmn.260909.getorg1770_rgj.mitos2110"))
        self.assertTrue(is_rescue_molecule("OG2998.ilmn.260909.getorg1770reseed_concat.mitos2110"))

    def test_real_molecules_are_not_rescue(self):
        self.assertFalse(is_rescue_molecule("getorg1770"))
        self.assertFalse(is_rescue_molecule("getorg1770reseed"))
        self.assertFalse(is_rescue_molecule("OG3000.ilmn.260909.getorg1770reseed.mitos2110"))
        self.assertFalse(is_rescue_molecule("OG2647.ilmn.260909.getorg1770.emma102"))
        self.assertFalse(is_rescue_molecule("mitohifi310"))

    def test_collapsed_is_a_real_molecule(self):
        """A collapse rewrites a confirmed concatemer to its monomer: that IS the molecule."""
        self.assertFalse(is_rescue_molecule("getorg1770_collapsed"))
        self.assertFalse(is_rescue_molecule("OG1.ilmn.260909.getorg1770_collapsed.mitos2110"))

    def test_empty_and_none(self):
        self.assertFalse(is_rescue_molecule(""))
        self.assertFalse(is_rescue_molecule(None))


class CrossLanguageDriftTests(unittest.TestCase):
    """The rule exists in Python and Groovy; a suffix added to one must reach the other.

    The Nextflow side cannot import the Python module, so the list is mirrored as
    isRescueMolecule() in subworkflows/local/ena_submission_prep/main.nf. Comparing the
    literals is the only thing that stops the two from quietly disagreeing -- and a
    disagreement here means an artificial molecule is withheld from webin but still reported
    submission_ready, or the reverse.
    """

    GROOVY = ROOT / "subworkflows" / "local" / "ena_submission_prep" / "main.nf"

    def test_groovy_and_python_suffix_lists_match(self):
        groovy = self.GROOVY.read_text()
        match = re.search(r"RESCUE_SUFFIXES\s*=\s*\[(.*?)\]", groovy, re.S)
        self.assertIsNotNone(match, f"RESCUE_SUFFIXES not found in {self.GROOVY}")
        groovy_suffixes = tuple(re.findall(r"'([^']+)'", match.group(1)))
        self.assertEqual(groovy_suffixes, RESCUE_SUFFIXES)

    def test_collapsed_is_absent_from_both(self):
        groovy = self.GROOVY.read_text()
        match = re.search(r"RESCUE_SUFFIXES\s*=\s*\[(.*?)\]", groovy, re.S)
        self.assertNotIn("_collapsed", match.group(1))
        self.assertNotIn("_collapsed", RESCUE_SUFFIXES)


if __name__ == "__main__":
    unittest.main()
