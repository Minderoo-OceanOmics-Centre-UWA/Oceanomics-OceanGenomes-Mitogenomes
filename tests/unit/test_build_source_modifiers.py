"""Unit tests for the lat_lon and Collection_date guards in build_source_modifiers.py.

An out-of-spec coordinate is not merely rejected downstream: EMBOSS demotes it to
a stray /note during EMBL conversion, so the bad value reaches ENA disguised as
free text rather than as an error. It has to be caught before the .src is written.

A bad collection_date fails harder and earlier -- table2asn raises
SEQ_DESCR.BadCollectionDate, which the validation gate counts as an ERROR and
quarantines the whole assembly for.

geo_loc_name resolution lives in bin/geo_loc_name_utils.py and is covered by
tests/unit/test_geo_loc_name_utils.py.
"""

import importlib.util
import sys
import types
import unittest
from datetime import date, timedelta
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]

# build_source_modifiers imports pandas and psycopg2 at module scope for the DB
# query path. The helper under test is pure, so stub the heavy imports rather than
# requiring a database client to run the unit suite. bin/ goes on sys.path so the
# geo_loc_name_utils sibling import resolves the same way it does under Nextflow,
# which bind-mounts bin/ into the task container and puts it on PATH.
sys.path.insert(0, str(ROOT / "bin"))
with mock.patch.dict(sys.modules, {
    "pandas": types.ModuleType("pandas"),
    "psycopg2": types.ModuleType("psycopg2"),
}):
    SPEC = importlib.util.spec_from_file_location(
        "build_source_modifiers", ROOT / "bin" / "build_source_modifiers.py"
    )
    MODULE = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(MODULE)


class ValidLatLonTests(unittest.TestCase):
    def test_well_formed_coordinate_is_accepted(self):
        self.assertTrue(MODULE.valid_lat_lon("18.14943 S 122.29572 E"))
        self.assertTrue(MODULE.valid_lat_lon("21.57628 S 156.49365 E"))

    def test_two_latitudes_are_rejected(self):
        # The real OG869 record: a hemisphere pair with no longitude at all.
        self.assertFalse(MODULE.valid_lat_lon("25.51000 N 23.43000 S"))

    def test_out_of_range_magnitudes_are_rejected(self):
        self.assertFalse(MODULE.valid_lat_lon("95.0 N 10.0 E"))
        self.assertFalse(MODULE.valid_lat_lon("18.1 S 200.0 E"))

    def test_empty_and_unknown_are_rejected(self):
        self.assertFalse(MODULE.valid_lat_lon(""))
        self.assertFalse(MODULE.valid_lat_lon("unknown"))

    def test_missing_hemisphere_is_rejected(self):
        self.assertFalse(MODULE.valid_lat_lon("18.14943 122.29572"))


class ValidCollectionDateTests(unittest.TestCase):
    def test_insdc_date_forms_are_accepted(self):
        self.assertTrue(MODULE.valid_collection_date("21-Oct-1952"))
        self.assertTrue(MODULE.valid_collection_date("Oct-1952"))
        self.assertTrue(MODULE.valid_collection_date("1952"))

    def test_single_digit_day_is_accepted(self):
        # pandas' %d zero-pads, but a hand-entered value may not.
        self.assertTrue(MODULE.valid_collection_date("3-Apr-2024"))

    def test_unknown_is_rejected(self):
        # The OG193 case: sample.date_collected was NULL and the old fillna wrote
        # this literal, which table2asn rejects as SEQ_DESCR.BadCollectionDate.
        self.assertFalse(MODULE.valid_collection_date("Unknown"))
        self.assertFalse(MODULE.valid_collection_date("unknown"))

    def test_insdc_missing_value_terms_are_rejected(self):
        # These are the tempting fix, and they do pass table2asn -- but they are
        # ENA *sample checklist* (ERC000011) vocabulary, not flatfile qualifier
        # values. Verified against ENA's own CollectionDateQualifierCheck in
        # sequencetools 2.33.2 (shipped in webin-cli 9.0.3): all three are
        # rejected, so using one would only move the failure to the last gate.
        self.assertFalse(MODULE.valid_collection_date("missing"))
        self.assertFalse(MODULE.valid_collection_date("not collected"))
        self.assertFalse(MODULE.valid_collection_date("not provided"))

    def test_future_date_is_rejected(self):
        # table2asn: "Collection_date is in the future"; ENA: FutureDateException.
        # The sample table holds day/month-transposed rows that land in the future.
        future = date.today() + timedelta(days=365)
        self.assertFalse(MODULE.valid_collection_date(future.strftime("%d-%b-%Y")))
        self.assertFalse(MODULE.valid_collection_date(str(future.year + 1)))

    def test_today_is_accepted(self):
        self.assertTrue(MODULE.valid_collection_date(date.today().strftime("%d-%b-%Y")))

    def test_empty_is_rejected(self):
        # Absent is handled by the caller's explicit emptiness test, not here.
        self.assertFalse(MODULE.valid_collection_date(""))
        self.assertFalse(MODULE.valid_collection_date("   "))

    def test_non_insdc_layouts_are_rejected(self):
        self.assertFalse(MODULE.valid_collection_date("2024-04-03"))
        self.assertFalse(MODULE.valid_collection_date("03/04/2024"))
        self.assertFalse(MODULE.valid_collection_date("21-Oct-52"))
        self.assertFalse(MODULE.valid_collection_date("32-Oct-1952"))
        self.assertFalse(MODULE.valid_collection_date("21-Xxx-1952"))


class FormatLatLonTests(unittest.TestCase):
    """The hemisphere is read, never invented.

    The sample table stores latitudes as unsigned magnitudes, so defaulting an
    unsigned value to the northern hemisphere placed fifteen Western Australian
    assemblies north of the equator. table2asn accepted the shape and rejected
    the value as SEQ_DESCR.LatLonValue, quarantining every one of them.
    """

    AU_NINGALOO = "Australia: Western Australia, Ningaloo"

    def test_unsigned_latitude_takes_its_hemisphere_from_the_country(self):
        # The exact batch-19 regression: OG2279, held as FAIL_TABLE2ASN.
        self.assertEqual(
            MODULE.format_lat_lon("22.03 113.891", self.AU_NINGALOO),
            ("22.03000 S 113.89100 E", "ok"),
        )

    def test_signed_latitude_wins_over_the_country(self):
        self.assertEqual(
            MODULE.format_lat_lon("-22.03 113.891", self.AU_NINGALOO),
            ("22.03000 S 113.89100 E", "ok"),
        )

    def test_stated_hemisphere_is_passed_through(self):
        self.assertEqual(
            MODULE.format_lat_lon("18.14943 S 122.29572 E", "Australia"),
            ("18.14943 S 122.29572 E", "ok"),
        )

    def test_northern_hemisphere_country_is_not_flipped(self):
        # The one locality in the table that is genuinely north of the equator.
        self.assertEqual(
            MODULE.format_lat_lon("29.5 34.9", "Israel: Elat, Gulf of Aquaba"),
            ("29.50000 N 34.90000 E", "ok"),
        )

    def test_dms_without_letters_uses_the_country(self):
        value, status = MODULE.format_lat_lon("22° 1.8 113° 53.5", self.AU_NINGALOO)
        self.assertEqual(status, "ok")
        self.assertTrue(value.startswith("22.03000 S 113.89"), value)

    def test_equator_straddling_country_cannot_place_a_bare_magnitude(self):
        self.assertEqual(
            MODULE.format_lat_lon("2.5 118.4", "Indonesia"), ("", "no_hemisphere")
        )

    def test_unmapped_country_cannot_place_a_bare_magnitude(self):
        self.assertEqual(
            MODULE.format_lat_lon("22.03 113.891", "Atlantis"), ("", "no_hemisphere")
        )

    def test_missing_value_country_cannot_place_a_bare_magnitude(self):
        self.assertEqual(
            MODULE.format_lat_lon("22.03 113.891", "not provided"),
            ("", "no_hemisphere"),
        )

    def test_stated_hemisphere_contradicting_the_country_is_dropped(self):
        # Exactly what table2asn raises SEQ_DESCR.LatLonValue for. Which of the
        # two is wrong cannot be told from here, so neither is corrected.
        self.assertEqual(
            MODULE.format_lat_lon("22.03 N 113.891 E", self.AU_NINGALOO),
            ("", "conflict"),
        )

    def test_two_latitudes_are_dropped(self):
        self.assertEqual(
            MODULE.format_lat_lon("25.51 N 23.43 S", "Australia"), ("", "unparsable")
        )

    def test_out_of_range_magnitude_is_dropped(self):
        self.assertEqual(
            MODULE.format_lat_lon("95.0 10.0", "Australia"), ("", "unparsable")
        )

    def test_absent_and_unknown_are_silent(self):
        for raw in (None, "", "   ", "unknown", "Unknown"):
            self.assertEqual(
                MODULE.format_lat_lon(raw, "Australia"), ("", "absent"), raw
            )

    def test_every_drop_status_has_an_operator_message(self):
        # A status with no entry in LATLON_DROP_REASONS would drop the coordinate
        # silently, which is the failure mode this whole path exists to avoid.
        for status in ("unparsable", "no_hemisphere", "conflict"):
            self.assertIn(status, MODULE.LATLON_DROP_REASONS)
        self.assertNotIn("absent", MODULE.LATLON_DROP_REASONS)
        self.assertNotIn("ok", MODULE.LATLON_DROP_REASONS)


class ParseCoordinateTests(unittest.TestCase):
    def test_bare_magnitude_reports_no_hemisphere(self):
        self.assertEqual(MODULE.parse_coordinate("22.03", "lat"), (22.03, None))

    def test_sign_is_a_recorded_hemisphere(self):
        self.assertEqual(MODULE.parse_coordinate("-22.03", "lat"), (22.03, "S"))
        self.assertEqual(MODULE.parse_coordinate("-113.891", "lon"), (113.891, "W"))

    def test_letter_is_read_and_magnitude_made_positive(self):
        self.assertEqual(MODULE.parse_coordinate("22.03 S", "lat"), (22.03, "S"))

    def test_letter_from_the_wrong_axis_is_rejected(self):
        # A longitude labelled S is not a longitude, and no fix is safe to guess.
        self.assertIsNone(MODULE.parse_coordinate("23.43 S", "lon"))
        self.assertIsNone(MODULE.parse_coordinate("23.43 E", "lat"))

    def test_non_numeric_is_rejected(self):
        self.assertIsNone(MODULE.parse_coordinate("not a coordinate", "lat"))


if __name__ == "__main__":
    unittest.main()
