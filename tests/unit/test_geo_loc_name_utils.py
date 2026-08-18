"""Unit tests for INSDC geo_loc_name resolution in geo_loc_name_utils.py.

geo_loc_name is one of only two mandatory fields in ENA checklist ERC000011, and
the sample table is rebuilt from a spreadsheet, so these mappings are the only
thing standing between a recorded typo and a rejected submission.
"""

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "geo_loc_name_utils", ROOT / "bin" / "geo_loc_name_utils.py"
)
geo_loc_name_utils = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(geo_loc_name_utils)

resolve = geo_loc_name_utils.resolve_geo_loc_name
CV = geo_loc_name_utils.INSDC_GEO_LOC_NAMES


class ControlledVocabularyTests(unittest.TestCase):
    def test_vocabulary_looks_like_the_insdc_list(self):
        self.assertEqual(len(CV), 294)
        for expected in ("Australia", "Palau", "Indian Ocean", "Southern Ocean"):
            self.assertIn(expected, CV)

    def test_every_alias_target_is_a_controlled_value(self):
        # An alias pointing at a non-CV value would swap one rejection for another.
        for target in geo_loc_name_utils.GEO_LOC_ALIASES.values():
            self.assertIn(target, CV, f"{target!r} is not in the INSDC vocabulary")

    def test_every_survey_target_is_a_controlled_value(self):
        for target in geo_loc_name_utils.SURVEY_GEO_LOC.values():
            self.assertIn(target, CV, f"{target!r} is not in the INSDC vocabulary")


class DirectAliasTests(unittest.TestCase):
    def test_long_form_political_name(self):
        self.assertEqual(
            resolve("Kingdom of Tonga: Tonga Trench"), ("Tonga: Tonga Trench", "aliased")
        )

    def test_long_form_name_without_locality(self):
        self.assertEqual(resolve("Republic of Palau"), ("Palau", "aliased"))

    def test_typos(self):
        self.assertEqual(resolve("Austalia: WA, Perth")[0], "Australia: WA, Perth")
        self.assertEqual(resolve("Kiribat: Kiritimati")[0], "Kiribati: Kiritimati")

    def test_casing_is_corrected(self):
        self.assertEqual(resolve("JAPAN: Okinawa"), ("Japan: Okinawa", "aliased"))

    def test_already_correct_value_is_untouched(self):
        self.assertEqual(resolve("Japan: Okinawa"), ("Japan: Okinawa", "ok"))
        self.assertEqual(
            resolve("Australia: WA, Thangoo"), ("Australia: WA, Thangoo", "ok")
        )

    def test_expanded_name_keeps_its_full_locality(self):
        self.assertEqual(
            resolve("Falklands: Stanley")[0],
            "Falkland Islands (Islas Malvinas): Stanley",
        )


class SurveyResolutionTests(unittest.TestCase):
    def test_high_seas_roo_rise_is_the_indian_ocean(self):
        self.assertEqual(
            resolve("High Seas: Roo Rise, Feature 1, South Wall"),
            ("Indian Ocean: Roo Rise, Feature 1, South Wall", "aliased"),
        )

    def test_chinese_research_vessel_is_the_pacific(self):
        value, status = resolve(
            "International Waters - Chinese Research Vessel: "
            "International Waters - Chinese Research Vessel, Portugal/China/Australia"
        )
        self.assertTrue(value.startswith("Pacific Ocean:"))
        self.assertEqual(status, "aliased")

    def test_south_shetland_is_the_southern_ocean(self):
        self.assertEqual(
            resolve("South Shetland: South Shetland"),
            ("Southern Ocean: South Shetland", "aliased"),
        )

    def test_same_label_at_a_different_survey_is_not_guessed(self):
        # 'High Seas' only means the Indian Ocean for the Roo Rise survey. Anywhere
        # else it has to be reported, not assumed.
        self.assertEqual(
            resolve("High Seas: Somewhere Else"),
            ("High Seas: Somewhere Else", "unmapped"),
        )


class MissingValueTests(unittest.TestCase):
    """geo_loc_name is mandatory, so 'no value' still has to be a legal value."""

    def test_nothing_recorded_becomes_the_insdc_missing_term(self):
        for value in ("Unknown", "unknown", "", "   ", None):
            self.assertEqual(resolve(value), ("not provided", "missing"))

    def test_missing_term_is_itself_a_controlled_value(self):
        # Substituting a term that is not in the vocabulary would swap a missing
        # mandatory field for a rejected one.
        self.assertIn(geo_loc_name_utils.MISSING_GEO_LOC, CV)

    def test_explicit_unknown_locality_is_also_missing(self):
        self.assertEqual(
            resolve("original locality unknown"), ("not provided", "missing")
        )


class BioSampleRecoveryTests(unittest.TestCase):
    """Countries recovered from BioSamples this project already registered."""

    SAMPLE_OG = "OG260"

    def test_map_is_populated_and_wholly_controlled(self):
        self.assertEqual(len(geo_loc_name_utils.BIOSAMPLE_GEO_LOC), 132)
        for og_id, value in geo_loc_name_utils.BIOSAMPLE_GEO_LOC.items():
            self.assertIn(value, CV, f"{og_id} maps to {value!r}, not an INSDC value")

    def test_missing_country_is_recovered_from_the_biosample(self):
        self.assertEqual(
            resolve("Unknown", self.SAMPLE_OG),
            (geo_loc_name_utils.BIOSAMPLE_GEO_LOC[self.SAMPLE_OG], "biosample"),
        )

    def test_recovery_applies_to_a_wholly_empty_value(self):
        self.assertEqual(resolve(None, self.SAMPLE_OG)[1], "biosample")

    def test_unknown_og_id_still_falls_back(self):
        self.assertEqual(resolve("Unknown", "OG999999"), ("not provided", "missing"))

    def test_recorded_country_always_wins_over_recovery(self):
        # The BioSample is only consulted when the row itself records nothing.
        self.assertEqual(
            resolve("Tonga: Tonga Trench", self.SAMPLE_OG), ("Tonga: Tonga Trench", "ok")
        )

    def test_recovery_never_rescues_an_unmapped_value(self):
        # An uncontrolled country must still quarantine so it gets corrected,
        # rather than being quietly replaced by the BioSample country.
        self.assertEqual(
            resolve("Narnia: Cair Paravel", self.SAMPLE_OG),
            ("Narnia: Cair Paravel", "unmapped"),
        )

    def test_og_id_is_optional(self):
        self.assertEqual(resolve("Unknown"), ("not provided", "missing"))


class DerivedFromLocalityTests(unittest.TestCase):
    """Rows with no country column arrive as their bare locality text."""

    def test_locality_naming_a_country_is_promoted(self):
        # The real OG184 record: no country, but the locality names one and its
        # coordinates (29 33 N, 034 57 E) agree.
        self.assertEqual(
            resolve("Israel, Elat, Gulf of Aquaba"),
            ("Israel: Elat, Gulf of Aquaba", "derived"),
        )

    def test_locality_naming_only_a_country(self):
        self.assertEqual(resolve("Israel"), ("Israel", "ok"))

    def test_locality_naming_nothing_controlled_is_not_guessed(self):
        self.assertEqual(
            resolve("Narnia, Cair Paravel"), ("Narnia, Cair Paravel", "unmapped")
        )

    def test_promotion_never_overrides_a_recorded_country(self):
        # A value that already has a country token is resolved on that token, so
        # the comma rule is never consulted.
        self.assertEqual(
            resolve("Australia: WA, Israel Bay"), ("Australia: WA, Israel Bay", "ok")
        )


class GuardTests(unittest.TestCase):

    def test_unmapped_value_is_returned_unchanged(self):
        # Returned as-is so table2asn flags it and the gate quarantines the sample.
        self.assertEqual(
            resolve("Narnia: Cair Paravel"), ("Narnia: Cair Paravel", "unmapped")
        )

    def test_resolution_is_idempotent(self):
        for value in (
            "Kingdom of Tonga: Tonga Trench",
            "High Seas: Roo Rise, Feature 2",
            "JAPAN: Okinawa",
            "Australia: WA, Perth",
            "Israel, Elat, Gulf of Aquaba",
        ):
            once, _ = resolve(value)
            self.assertEqual(resolve(once), (once, "ok"))

    def test_every_resolved_value_starts_with_a_controlled_token(self):
        # The whole point: nothing that resolves may still be rejected by ENA.
        for value in (
            "Australia: WA, Thangoo, Roebuck Bay",
            "High Seas: Roo Rise, Feature 1, South Wall",
            "Republic of Palau: Ngerukewid",
            "JAPAN: Okinawa",
            "Austalia: WA, Perth",
            "Kiribat: Kiritimati",
            "Falklands: Stanley",
            "Kingdom of Tonga: Tonga Trench",
            "South Shetland: South Shetland",
        ):
            resolved, status = resolve(value)
            self.assertNotEqual(status, "unmapped")
            self.assertIn(resolved.partition(":")[0].strip(), CV)


class UnmappedWarningTests(unittest.TestCase):
    def test_warning_names_the_line_to_add(self):
        message = geo_loc_name_utils.unmapped_warning("OG1234", "Repulic of Palau: X")
        self.assertIn("OG1234", message)
        self.assertIn("'Repulic of Palau'", message)
        self.assertIn('"repulic of palau":', message)


if __name__ == "__main__":
    unittest.main()
