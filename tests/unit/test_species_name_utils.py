"""Unit tests for open-nomenclature normalisation in species_name_utils.py."""

import importlib.util
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "species_name_utils", ROOT / "bin" / "species_name_utils.py"
)
species_name_utils = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(species_name_utils)

normalise = species_name_utils.normalise_open_nomenclature


class OpenNomenclatureTests(unittest.TestCase):
    def test_missing_period_is_added(self):
        # The OG1834 failure: ENA rejects "Chaunax sp" as not submittable.
        self.assertEqual(normalise("Chaunax sp"), "Chaunax sp.")
        self.assertEqual(normalise("Hoplichthys sp"), "Hoplichthys sp.")

    def test_spp_collapses_to_sp(self):
        self.assertEqual(normalise("Blachea spp."), "Blachea sp.")
        self.assertEqual(normalise("Astronesthes spp."), "Astronesthes sp.")
        self.assertEqual(normalise("Zenion spp"), "Zenion sp.")

    def test_idempotent(self):
        self.assertEqual(normalise("Chaunax sp."), "Chaunax sp.")
        self.assertEqual(normalise(normalise("Chaunax sp")), "Chaunax sp.")

    def test_binomials_untouched(self):
        for name in (
            "Lutjanus quinquelineatus",
            "Limnichthys fasciatus",
            "Scobinichthys granulatus",
            "Chaunax penicillatus",
        ):
            self.assertEqual(normalise(name), name)

    def test_out_of_scope_forms_pass_through(self):
        # Deliberately not handled: these need a judgement call, not a regex.
        for name in (
            "Centrodraco sp 2",
            "Hymenocephalus sp 4",
            "Nesogobius sp. `groove cheek`",
            "Eptatretus sp (Eptatretus cf. goliath)",
            "Synodus macrops cf",
            "Neobythites bimarginatus cf",
        ):
            self.assertEqual(normalise(name), name)

    def test_surrounding_whitespace_is_stripped(self):
        self.assertEqual(normalise("  Chaunax sp  "), "Chaunax sp.")

    def test_empty_and_unknown_pass_through(self):
        self.assertEqual(normalise(""), "")
        self.assertIsNone(normalise(None))
        self.assertEqual(normalise("unknown"), "unknown")


class InlineModuleCopyTests(unittest.TestCase):
    """
    modules/local/validated_species_query/main.nf can't import the helper (its
    script is a heredoc running from the task work dir), so it inlines the same
    rule. Guard against the two drifting apart.
    """

    MODULE = ROOT / "modules" / "local" / "validated_species_query" / "main.nf"

    def test_inline_regex_matches_helper_behaviour(self):
        text = self.MODULE.read_text()
        found = re.search(r're\.fullmatch\(r"(.+?)"', text)
        self.assertIsNotNone(found, "inline normalisation regex not found in module")
        # The .nf is a Groovy string, so backslashes are doubled in the source.
        inline = re.compile(found.group(1).replace("\\\\", "\\"))

        for name in ("Chaunax sp", "Blachea spp.", "Zenion spp", "Chaunax sp."):
            match = inline.fullmatch(name.strip())
            self.assertIsNotNone(match, f"inline regex failed to match {name!r}")
            self.assertEqual(match.group(1) + " sp.", normalise(name))

        for name in (
            "Lutjanus quinquelineatus",
            "Centrodraco sp 2",
            "Synodus macrops cf",
            "unknown",
        ):
            self.assertIsNone(
                inline.fullmatch(name.strip()),
                f"inline regex should not match {name!r}",
            )


if __name__ == "__main__":
    unittest.main()
