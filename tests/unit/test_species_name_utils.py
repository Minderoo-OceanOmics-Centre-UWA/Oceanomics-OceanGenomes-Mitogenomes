"""Unit tests for open-nomenclature normalisation in species_name_utils.py."""

import importlib.util
import sys
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
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
    modules/local/validated_species_query/main.nf used to INLINE a copy of the
    normalisation rule, on the grounds that its heredoc runs from the task work
    dir and so cannot do a sibling import. That was true, but the copy then
    drifted in scope: it handled only 'sp'/'spp' and left a BARE GENUS alone, so
    the qc-only path still emitted an organism ENA rejects as not submittable
    after the sample had cleared every gate upstream.

    It now locates bin/ on PATH and imports the canonical helper. These tests pin
    that: the copy must stay gone, and the import must stay.
    """

    MODULE = ROOT / "modules" / "local" / "validated_species_query" / "main.nf"

    def test_the_module_imports_the_canonical_helper(self):
        text = self.MODULE.read_text()
        self.assertIn("from species_name_utils import normalise_open_nomenclature", text)
        # And finds bin/ the only way it can from a task work dir.
        self.assertIn('os.environ.get("PATH"', text)

    def test_no_inline_copy_of_the_normalisation_rule_remains(self):
        text = self.MODULE.read_text()
        self.assertIsNone(
            re.search(r're\.fullmatch\(r"\(\[A-Za-z\]', text),
            "an inline copy of the normalisation regex has come back; import the "
            "helper instead -- two implementations of one rule is what drifted before",
        )

    def test_the_helper_covers_what_the_old_inline_copy_missed(self):
        # The specific gap that made the copy a bug rather than a duplication.
        self.assertEqual(normalise("Serrivomer"), "Serrivomer sp.")


parse_nominal = species_name_utils.parse_nominal


class BareGenusNormalisationTests(unittest.TestCase):
    """A bare genus is not a submittable organism, and it used to reach ENA.

    'Serrivomer' passed the QC gate only because the old BLAST test was a
    substring search that matched inside 'Serrivomer jesperseni', and was then
    rejected at Webin: ENA recognises 'Serrivomer sp.' and does not recognise
    'Serrivomer'. This half of the fix changes no gate decision at all.
    """

    def test_a_bare_genus_becomes_genus_sp(self):
        self.assertEqual(normalise("Serrivomer"), "Serrivomer sp.")
        self.assertEqual(normalise("  Zenion  "), "Zenion sp.")

    def test_it_is_idempotent(self):
        self.assertEqual(normalise(normalise("Serrivomer")), "Serrivomer sp.")

    def test_a_family_name_is_not_turned_into_a_species(self):
        # That decision belongs to the family-rank handling in the validator, not
        # to the normaliser: 'Ophidiidae sp.' asserts an undescribed species in a
        # genus called Ophidiidae, which is not what the label says.
        for name in ("Ophidiidae", "Macrouridae", "Serraninae"):
            self.assertEqual(normalise(name), name)

    def test_unresolved_sentinels_pass_through_untouched(self):
        # 'unknown sp.' is worse than 'unknown', because it looks submittable.
        for name in ("unknown", "NA", "n/a", "none", "TBC", ""):
            self.assertEqual(normalise(name), name)

    def test_real_binomials_are_untouched(self):
        for name in ("Notacanthus abbotti", "Chrysophrys auratus"):
            self.assertEqual(normalise(name), name)


class ParseNominalTests(unittest.TestCase):
    """Classify a nominal ID by the rank it actually asserts."""

    def test_a_binomial_asserts_a_species(self):
        self.assertEqual(parse_nominal("Notacanthus abbotti"),
                         ("species", "Notacanthus abbotti"))

    def test_open_nomenclature_asserts_a_genus(self):
        for name in ("Diaphus sp.", "Diaphus sp", "Diaphus spp", "Diaphus spp.",
                     "Diaphus sp 1", "Diaphus sp. 2"):
            self.assertEqual(parse_nominal(name), ("genus", "Diaphus"), name)

    def test_a_bare_genus_asserts_a_genus(self):
        self.assertEqual(parse_nominal("Serrivomer"), ("genus", "Serrivomer"))

    def test_a_family_name_asserts_a_family(self):
        self.assertEqual(parse_nominal("Ophidiidae"), ("family", "Ophidiidae"))
        self.assertEqual(parse_nominal("Ophidiidae sp."), ("family", "Ophidiidae"))

    def test_a_family_name_wins_over_a_trailing_field_note(self):
        # 'Macrouridae black' is a family and a colour, not a species. Reading it
        # as a binomial is what sent a reference lookup searching for a species
        # that does not exist, and took that sample out of the run entirely.
        self.assertEqual(parse_nominal("Macrouridae black"),
                         ("family", "Macrouridae"))

    def test_uncertainty_markers_assert_a_tentative_species(self):
        for name in ("Squalus notocaudatus?", "Synodus macrops cf",
                     "Synodus cf. macrops", "Squalus aff. megalops"):
            rank, value = parse_nominal(name)
            self.assertEqual(rank, "species_uncertain", name)
            self.assertIn(value, ("Squalus", "Synodus"), name)

    def test_parenthetical_field_notes_are_stripped(self):
        self.assertEqual(parse_nominal("Diaphus sp 2 (short jaw group)"),
                         ("genus", "Diaphus"))
        self.assertEqual(parse_nominal("Spectrunculus grandis (Brown color)"),
                         ("species", "Spectrunculus grandis"))

    def test_unresolved_values_assert_nothing(self):
        for name in (None, "", "   ", "unknown", "NA", "n/a", "none", "TBC"):
            self.assertEqual(parse_nominal(name), (None, None), repr(name))

    def test_an_order_or_class_name_under_fires_to_genus_which_is_safe(self):
        # Documented limitation: the family heuristic keys on -idae/-inae, so an
        # order or class name (which invertebrate labels carry far more often than
        # fish ones) classifies as a genus. That can only FAIL to match -- no BLAST
        # genus set or LCA genus column carries 'Scleractinia' -- so the sample
        # stays held, which is the safe direction.
        self.assertEqual(parse_nominal("Scleractinia"), ("genus", "Scleractinia"))
