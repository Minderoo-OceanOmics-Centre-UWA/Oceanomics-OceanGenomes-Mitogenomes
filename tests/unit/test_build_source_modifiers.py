"""Unit tests for the lat_lon guard in build_source_modifiers.py.

An out-of-spec coordinate is not merely rejected downstream: EMBOSS demotes it to
a stray /note during EMBL conversion, so the bad value reaches ENA disguised as
free text rather than as an error. It has to be caught before the .src is written.

geo_loc_name resolution lives in bin/geo_loc_name_utils.py and is covered by
tests/unit/test_geo_loc_name_utils.py.
"""

import importlib.util
import sys
import types
import unittest
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


if __name__ == "__main__":
    unittest.main()
