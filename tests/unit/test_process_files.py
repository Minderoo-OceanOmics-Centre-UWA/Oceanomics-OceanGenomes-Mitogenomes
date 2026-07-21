"""Unit tests for filename-derived metadata in process_files.py."""

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "process_files", ROOT / "bin" / "process_files.py"
)
process_files = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(process_files)


class AssemblyMethodTests(unittest.TestCase):
    def test_getorganelle(self):
        self.assertEqual(
            process_files.derive_assembly_method("OG1.ilmn.260101.v177getorg.emma102"),
            "GetOrganelle v.1.7.7",
        )

    def test_mitohifi(self):
        self.assertEqual(
            process_files.derive_assembly_method("OG1.hifi.260101.v323mitohifi.emma102"),
            "MitoHifi v.3.2.3",
        )

    def test_oatk(self):
        self.assertEqual(
            process_files.derive_assembly_method("OG1422.hifi.260227.v10oatk.emma102"),
            "Oatk v.1.0",
        )

    def test_unknown_assembler_is_rejected(self):
        with self.assertRaisesRegex(SystemExit, "Unknown assembler code"):
            process_files.derive_assembly_method("OG1.hifi.260101.v10unknown.emma102")


if __name__ == "__main__":
    unittest.main()
