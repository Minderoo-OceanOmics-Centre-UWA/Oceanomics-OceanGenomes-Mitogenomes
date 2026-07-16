import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "compile_upload_report", ROOT / "bin" / "compile_upload_report.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class EnaUploadReportTests(unittest.TestCase):
    def test_ena_status_classification(self):
        self.assertEqual(
            MODULE.classify_status("ena_validation", "✅ Success: inserted ENA validation attempt\nUPLOAD_EXIT=0\n"),
            "success",
        )
        self.assertEqual(
            MODULE.classify_status("ena_validation", "⚠️ Exact ENA validation result already recorded\nUPLOAD_EXIT=0\n"),
            "preserved",
        )
        self.assertEqual(
            MODULE.classify_status("ena_validation", "❌ Database error: unavailable\nUPLOAD_EXIT=1\n"),
            "failed",
        )

    def test_assembly_prefix_parsing_and_summary_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            upload = root / "OG1.hifi.260101.final.ena_validation.upload.txt"
            upload.write_text("✅ Success: inserted ENA validation attempt\nUPLOAD_EXIT=0\n")
            prefix, og_id, step = MODULE.parse_filename(upload)
            self.assertEqual(prefix, "OG1.hifi.260101.final")
            self.assertEqual(og_id, "OG1")
            self.assertEqual(step, "ena_validation")
            output = root / "summary.tsv"
            MODULE.write_summary_tsv(MODULE.collect_inputs(root), output)
            lines = output.read_text().splitlines()
            self.assertIn("ena_validation", lines[0].split("\t"))
            self.assertEqual(lines[1].split("\t")[-1], "success")


if __name__ == "__main__":
    unittest.main()
