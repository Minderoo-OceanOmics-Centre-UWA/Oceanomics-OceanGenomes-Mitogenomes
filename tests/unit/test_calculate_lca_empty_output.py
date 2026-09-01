"""A region with no valid BLAST hits must still produce header-only LCA files.

Downstream, upload_results_mito closes each sample's result group once its
region count is reached, so a hit-less region that silently produced no file
would leave the group one short and stall the sample until the end of the run.
"""

import ast
import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
BIN = ROOT / "bin"
ASSETS = ROOT / "assets"


def load_calculate_lca():
    """Import calculateLCA.py without requiring pandas to be installed."""
    spec = importlib.util.spec_from_file_location(
        "calculateLCA_under_test", BIN / "calculateLCA.py"
    )
    module = importlib.util.module_from_spec(spec)
    try:
        import pandas  # noqa: F401

        patches = {}
    except ImportError:
        # calculateLCA.py uses pd.DataFrame in type annotations evaluated at
        # class-definition time, so the stand-in needs that attribute to exist.
        stub = types.ModuleType("pandas")
        stub.DataFrame = object
        stub.Series = object
        patches = {"pandas": stub}
    with mock.patch.dict(sys.modules, patches):
        spec.loader.exec_module(module)
    return module


class WriteResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_calculate_lca()

    def _writer(self):
        # write_results only touches self.logger, so a bare instance is enough.
        analyser = object.__new__(self.module.BLASTLCAAnalyzer)
        analyser.logger = mock.MagicMock()
        return analyser

    def test_empty_results_write_header_only_file(self):
        writer = self._writer()
        for columns in (self.module.TAXA_RAW_COLUMNS, self.module.TAXA_FINAL_COLUMNS):
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "out.tsv"
                writer.write_results([], out, columns)
                lines = out.read_text().splitlines()
                self.assertEqual(len(lines), 1)
                self.assertEqual(lines[0].split("\t"), list(columns))

    def test_populated_results_keep_their_own_header(self):
        writer = self._writer()
        rows = [{"seq_id": "a", "species_in_LCA": "Genus species"}]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.tsv"
            writer.write_results(rows, out, self.module.TAXA_FINAL_COLUMNS)
            lines = out.read_text().splitlines()
            self.assertEqual(lines[0], "seq_id\tspecies_in_LCA")
            self.assertEqual(lines[1], "a\tGenus species")

    def test_no_columns_and_no_rows_is_an_error(self):
        writer = self._writer()
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                writer.write_results([], Path(tmp) / "out.tsv")

    def test_species_in_lca_column_is_present_for_the_empty_case(self):
        # species_validation.py reads this column by name from the combined LCA
        # file; a header missing it would make every zero-hit sample unparseable.
        self.assertIn("species_in_LCA", self.module.TAXA_FINAL_COLUMNS)


class ColumnConstantDriftTests(unittest.TestCase):
    """The header constants must stay in step with the rows actually written.

    TAXA_*_COLUMNS only supply the header for the empty case, so a new field
    added to the taxaRaw / taxaFinal dicts would not fail any normal run -- it
    would just silently give hit-less regions a header the populated files no
    longer share. Compare the two directly instead of trusting review.
    """

    @classmethod
    def setUpClass(cls):
        cls.tree = ast.parse((BIN / "calculateLCA.py").read_text())

    def _constant(self, name):
        for node in self.tree.body:
            if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                if node.targets[0].id == name:
                    return [element.value for element in node.value.elts]
        self.fail(f"{name} not found in calculateLCA.py")

    def _appended_dict_keys(self, list_name):
        for node in ast.walk(self.tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "append"
                and getattr(node.func.value, "id", None) == list_name
                and node.args
                and isinstance(node.args[0], ast.Dict)
            ):
                return [key.value for key in node.args[0].keys]
        self.fail(f"{list_name}.append({{...}}) not found in calculateLCA.py")

    def test_raw_columns_match_the_written_rows(self):
        self.assertEqual(
            self._constant("TAXA_RAW_COLUMNS"), self._appended_dict_keys("taxaRaw")
        )

    def test_final_columns_match_the_written_rows(self):
        self.assertEqual(
            self._constant("TAXA_FINAL_COLUMNS"), self._appended_dict_keys("taxaFinal")
        )


class EmptyLcaAssetTests(unittest.TestCase):
    """The placeholder asset used for samples that annotated zero regions."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_calculate_lca()

    def test_asset_header_matches_the_real_lca_header(self):
        header = (ASSETS / "empty_lca.tsv").read_text().splitlines()
        self.assertEqual(len(header), 1, "placeholder must be header-only")
        self.assertEqual(header[0].split("\t"), list(self.module.TAXA_FINAL_COLUMNS))

    def test_empty_blast_placeholder_is_empty(self):
        # Filtered BLAST output is headerless, so any content here would be
        # parsed as a hit row by species_validation.load_blast_species_set.
        self.assertEqual((ASSETS / "empty_blast_filtered.tsv").read_text(), "")


if __name__ == "__main__":
    unittest.main()
