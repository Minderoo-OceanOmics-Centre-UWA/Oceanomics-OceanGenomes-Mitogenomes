"""Tests for the dependency-light mean_depth reader.

The value this returns becomes COVERAGE in an ENA manifest, so the distinction
that matters most is None (unmeasured, COVERAGE omitted) versus 0.0 (measured as
zero, COVERAGE emitted). Conflating them states something false about the
assembly rather than staying silent about it.

Kept in step with parse_depth_tsv in push_mtdna_assm_results.py, which is the
pandas-based canonical reader; the two must agree on what counts as unusable.
"""

import importlib.util
import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("depth_tsv", ROOT / "bin" / "depth_tsv.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

PLACEHOLDER = ROOT / "assets" / "placeholders" / "empty_mito_depth.tsv"


class ReadMeanDepthTests(unittest.TestCase):
    def _tsv(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write(text)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_a_real_measurement_is_returned_as_a_float(self):
        self.assertEqual(
            MODULE.read_mean_depth(self._tsv("sample\tmean_depth\nOG1\t812.5\n")), 812.5
        )

    def test_a_measured_zero_survives_as_zero(self):
        """Must not be swallowed as falsy: zero coverage is a measurement."""
        value = MODULE.read_mean_depth(self._tsv("sample\tmean_depth\nOG1\t0\n"))
        self.assertIsNotNone(value)
        self.assertEqual(value, 0.0)

    def test_the_real_placeholder_reads_as_unmeasured(self):
        """The exact file the pipeline substitutes when there is no depth."""
        self.assertIsNone(MODULE.read_mean_depth(str(PLACEHOLDER)))

    def test_no_path_no_file_and_no_columns_are_all_unmeasured(self):
        self.assertIsNone(MODULE.read_mean_depth(None))
        self.assertIsNone(MODULE.read_mean_depth("/nonexistent/mito_depth.tsv"))
        self.assertIsNone(MODULE.read_mean_depth(self._tsv("")))
        self.assertIsNone(MODULE.read_mean_depth(self._tsv("sample\tbreadth_1x\nOG1\t1\n")))

    def test_blank_na_and_unparseable_cells_are_unmeasured(self):
        for cell in ("", "  ", "NA", "na", "n/a", "unknown"):
            with self.subTest(cell=cell):
                self.assertIsNone(
                    MODULE.read_mean_depth(self._tsv(f"sample\tmean_depth\nOG1\t{cell}\n"))
                )

    def test_only_the_first_data_row_is_read(self):
        """mito_depth.py writes one row per assembly; matches df.iloc[0]."""
        self.assertEqual(
            MODULE.read_mean_depth(
                self._tsv("sample\tmean_depth\nOG1\t100\nOG2\t900\n")
            ),
            100.0,
        )


if __name__ == "__main__":
    unittest.main()
