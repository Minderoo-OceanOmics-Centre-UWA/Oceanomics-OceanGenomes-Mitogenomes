"""Tests for persistent, atomic LCA cache preparation."""

import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

try:
    import pandas as pd
    import pyarrow  # noqa: F401
except ImportError:
    pd = None


ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"


class CacheSignatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(BIN))
        spec = importlib.util.spec_from_file_location(
            "prepare_lca_databases_signature", BIN / "prepare_lca_databases.py"
        )
        cls.module = importlib.util.module_from_spec(spec)
        fake_calculate_lca = types.ModuleType("calculateLCA")
        fake_calculate_lca.Config = object
        fake_calculate_lca.DatabaseManager = object
        fake_calculate_lca.FISHBASE_DATABASES = {}
        with mock.patch.dict(
            sys.modules,
            {"pandas": types.ModuleType("pandas"), "calculateLCA": fake_calculate_lca},
        ):
            spec.loader.exec_module(cls.module)

    @classmethod
    def tearDownClass(cls):
        if sys.path and sys.path[0] == str(BIN):
            sys.path.pop(0)

    def test_signature_is_order_independent(self):
        first = {
            "b.parquet": {"bytes": 2, "sha256": "b" * 64},
            "a.parquet": {"bytes": 1, "sha256": "a" * 64},
        }
        second = dict(reversed(list(first.items())))
        self.assertEqual(
            self.module.cache_signature(first),
            self.module.cache_signature(second),
        )

    def test_signature_includes_schema_and_file_metadata(self):
        files = {"a.parquet": {"bytes": 1, "sha256": "a" * 64}}
        initial = self.module.cache_signature(files)
        changed_file = self.module.cache_signature(
            {"a.parquet": {"bytes": 2, "sha256": "b" * 64}}
        )
        with mock.patch.object(self.module, "CACHE_SCHEMA_VERSION", 2):
            changed_schema = self.module.cache_signature(files)
        self.assertNotEqual(initial, changed_file)
        self.assertNotEqual(initial, changed_schema)


@unittest.skipUnless(pd is not None, "pandas/pyarrow not available")
class PrepareLCADatabasesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(BIN))
        spec = importlib.util.spec_from_file_location(
            "prepare_lca_databases", BIN / "prepare_lca_databases.py"
        )
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
        import calculateLCA
        cls.calculate_lca = calculateLCA

    @classmethod
    def tearDownClass(cls):
        if sys.path and sys.path[0] == str(BIN):
            sys.path.pop(0)

    def make_valid_cache(self, root):
        pd.DataFrame({
            "SpecCode": [1], "Genus": ["Testus"], "Species": ["fish"], "FamCode": [1]
        }).to_parquet(root / "fishbase_species.parquet")
        pd.DataFrame({
            "FamCode": [1], "Family": ["Testidae"], "Order": ["Testiformes"], "Class": ["Actinopteri"]
        }).to_parquet(root / "fishbase_families.parquet")
        pd.DataFrame({
            "SpecCode": [1], "SynGenus": ["Oldtestus"], "SynSpecies": ["fish"]
        }).to_parquet(root / "fishbase_synonyms.parquet")
        taxdump = root / "taxdump"
        taxdump.mkdir()
        (taxdump / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
        (taxdump / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")

    def run_prepare(self, cache, output_manifest):
        with mock.patch.object(self.module.DatabaseManager, "load_fishbase_data"), \
             mock.patch.object(self.module.DatabaseManager, "load_ncbi_taxdump"):
            self.module.prepare(cache, output_manifest)
        return json.loads(output_manifest.read_text())

    def test_prepare_emits_verified_checksum_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache"
            cache.mkdir()
            self.make_valid_cache(cache)
            output_manifest = root / "task_manifest.json"

            task_manifest = self.run_prepare(cache, output_manifest)
            cache_manifest = json.loads((cache / "manifest.json").read_text())
            self.assertEqual(task_manifest, cache_manifest)
            self.assertEqual(task_manifest["status"], "verified")
            self.assertRegex(task_manifest["cache_signature"], r"^[0-9a-f]{64}$")
            self.assertEqual(len(task_manifest["files"]), 5)
            self.assertTrue(all(entry["sha256"] for entry in task_manifest["files"].values()))

    def test_signature_is_stable_across_verification_times(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache"
            cache.mkdir()
            self.make_valid_cache(cache)

            first = self.run_prepare(cache, root / "first.json")
            second = self.run_prepare(cache, root / "second.json")

            self.assertNotEqual(first["verified_at"], second["verified_at"])
            self.assertEqual(first["cache_signature"], second["cache_signature"])
            self.assertEqual(second, json.loads((cache / "manifest.json").read_text()))

    def test_signature_changes_with_file_or_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache"
            cache.mkdir()
            self.make_valid_cache(cache)

            initial = self.run_prepare(cache, root / "initial.json")
            names = cache / "taxdump" / "names.dmp"
            names.write_text(names.read_text() + "2\t|\ttest\t|\t\t|\tscientific name\t|\n")
            changed_file = self.run_prepare(cache, root / "changed-file.json")
            self.assertNotEqual(initial["cache_signature"], changed_file["cache_signature"])

            with mock.patch.object(self.module, "CACHE_SCHEMA_VERSION", 2):
                changed_schema = self.run_prepare(cache, root / "changed-schema.json")
            self.assertNotEqual(changed_file["cache_signature"], changed_schema["cache_signature"])

    def test_invalid_cache_files_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)

            for case in ("missing", "empty", "invalid-columns"):
                with self.subTest(case=case):
                    cache = root / case
                    cache.mkdir()
                    self.make_valid_cache(cache)
                    species = cache / "fishbase_species.parquet"
                    if case == "missing":
                        species.unlink()
                    elif case == "empty":
                        species.write_bytes(b"")
                    else:
                        pd.DataFrame({"wrong": [1]}).to_parquet(species)

                    with self.assertRaises((RuntimeError, OSError)):
                        self.run_prepare(cache, root / f"{case}.json")

    def test_fishbase_download_retries_and_writes_atomically(self):
        with tempfile.TemporaryDirectory() as tmp:
            manager = self.calculate_lca.DatabaseManager(Path(tmp))
            expected = pd.DataFrame({"value": [1]})
            with mock.patch.object(
                self.calculate_lca.pd,
                "read_parquet",
                side_effect=[OSError("503"), OSError("504"), expected, expected],
            ) as reader, mock.patch.object(self.calculate_lca.time, "sleep"):
                observed = manager._download_with_cache("https://example.invalid/data.parquet", "data.parquet")

            self.assertEqual(reader.call_count, 4)
            self.assertEqual(observed.to_dict(), expected.to_dict())
            self.assertTrue((Path(tmp) / "data.parquet").is_file())
            self.assertEqual(list(Path(tmp).glob("*.part")), [])


if __name__ == "__main__":
    unittest.main()
