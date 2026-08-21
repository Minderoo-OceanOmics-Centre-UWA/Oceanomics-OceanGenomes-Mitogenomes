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
        self.assertEqual(
            MODULE.classify_status(
                "ena_validation",
                "🔁 Updated ENA validation attempt for OG1.ilmn.240313.getorg1770.emma102 "
                "(overwrote the previous attempt)\nUPLOAD_EXIT=0\n",
            ),
            "success_updated",
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
            header, row = (line.split("\t") for line in lines)
            self.assertIn("ena_validation", header)
            self.assertEqual(dict(zip(header, row))["ena_validation"], "success")

    def _write(self, root, name, text):
        (root / name).write_text(text)

    def test_ena_annotation_suffix_folds_onto_the_assembly_row(self):
        """The ENA file is named after full_seqid, one token longer than its siblings."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root, "OG1.hic.250624.v3mitohifi.mtdna.upload.txt",
                        "\u2705 Success: Inserted/Updated mitogenome_data\n")
            self._write(root, "OG1.hic.250624.v3mitohifi.emma102.ena_validation.upload.txt",
                        "\u2705 Success: inserted ENA validation attempt\nUPLOAD_EXIT=0\n")
            grouped = MODULE.collect_inputs(root)
            self.assertEqual(list(grouped), ["OG1.hic.250624.v3mitohifi"])
            entry = grouped["OG1.hic.250624.v3mitohifi"]
            self.assertEqual(entry["annotation"], "emma102")
            self.assertEqual(set(entry["steps"]), {"assembly", "ena_validation"})

            output = root / "summary.tsv"
            MODULE.write_summary_tsv(grouped, output)
            header, row = (line.split("\t") for line in output.read_text().splitlines())
            values = dict(zip(header, row))
            self.assertEqual(values["annotation_version"], "emma102")
            self.assertEqual(values["assembly"], "success")
            self.assertEqual(values["ena_validation"], "success")

    def test_ena_only_input_falls_back_to_the_four_field_prefix(self):
        """The standalone ena.nf entry point produces no anchor to match against."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root, "OG1.ilmn.240313.getorg1770.emma102.ena_validation.upload.txt",
                        "\u2705 Success: inserted ENA validation attempt\nUPLOAD_EXIT=0\n")
            grouped = MODULE.collect_inputs(root)
            self.assertEqual(list(grouped), ["OG1.ilmn.240313.getorg1770"])
            self.assertEqual(grouped["OG1.ilmn.240313.getorg1770"]["annotation"], "emma102")

    def test_two_annotation_versions_keep_the_worse_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root, "OG1.hic.250624.v3mitohifi.mtdna.upload.txt",
                        "\u2705 Success: Inserted/Updated mitogenome_data\n")
            self._write(root, "OG1.hic.250624.v3mitohifi.emma102.ena_validation.upload.txt",
                        "\u2705 Success: inserted ENA validation attempt\nUPLOAD_EXIT=0\n")
            self._write(root, "OG1.hic.250624.v3mitohifi.emma103.ena_validation.upload.txt",
                        "\u274c Database error: unavailable\nUPLOAD_EXIT=1\n")
            grouped = MODULE.collect_inputs(root)
            self.assertEqual(list(grouped), ["OG1.hic.250624.v3mitohifi"])
            entry = grouped["OG1.hic.250624.v3mitohifi"]
            self.assertEqual(entry["annotation"], "emma102,emma103")
            output = root / "summary.tsv"
            MODULE.write_summary_tsv(grouped, output)
            header, row = (line.split("\t") for line in output.read_text().splitlines())
            self.assertEqual(dict(zip(header, row))["ena_validation"], "failed")


class QcValidatorReportTests(unittest.TestCase):
    def _write(self, root, name, text):
        (root / name).write_text(text)

    def test_qc_validator_status_classification(self):
        self.assertEqual(
            MODULE.classify_status(
                "qc_validator",
                "\u2705 Success: lca_validation validator_2 set to 'QCd-nf-core' "
                "for OG1.hifi.260101.v3mitohifi.emma102\nUPLOAD_EXIT=0\n",
            ),
            "validator_2_set",
        )
        self.assertEqual(
            MODULE.classify_status(
                "qc_validator",
                "\u26a0\ufe0f Existing validator_2 preserved for "
                "OG1.hifi.260101.v3mitohifi.emma102: validator_2='TP'\nUPLOAD_EXIT=0\n",
            ),
            "preserved",
        )
        self.assertEqual(
            MODULE.classify_status(
                "qc_validator",
                "\u2139\ufe0f OG1.hifi.260101.v3mitohifi.emma102 not submission_ready "
                "\u2014 skipping validator_2 write.\nUPLOAD_EXIT=0\n",
            ),
            "not_ready",
        )
        self.assertEqual(
            MODULE.classify_status(
                "qc_validator",
                "\u26a0\ufe0f No lca_validation row for OG1.hifi.260101.v3mitohifi.emma102 "
                "\u2014 validator_2 not set.\nUPLOAD_EXIT=0\n",
            ),
            "no_row",
        )
        self.assertEqual(
            MODULE.classify_status(
                "qc_validator",
                "\u274c Database error: unavailable\nUPLOAD_EXIT=1\n",
            ),
            "failed",
        )

    def test_qc_validator_shares_the_assembly_row_with_ena_validation(self):
        """Both files are named after full_seqid, so both need the annotation stripped."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write(root, "OG1.hic.250624.v3mitohifi.mtdna.upload.txt",
                        "\u2705 Success: Inserted/Updated mitogenome_data\n")
            self._write(root, "OG1.hic.250624.v3mitohifi.emma102.ena_validation.upload.txt",
                        "\u2705 Success: inserted ENA validation attempt\nUPLOAD_EXIT=0\n")
            self._write(root, "OG1.hic.250624.v3mitohifi.emma102.qc_validator.upload.txt",
                        "\u2705 Success: lca_validation validator_2 set to 'QCd-nf-core' "
                        "for OG1.hic.250624.v3mitohifi.emma102\nUPLOAD_EXIT=0\n")
            grouped = MODULE.collect_inputs(root)
            # One row, not one per full_seqid-keyed step.
            self.assertEqual(list(grouped), ["OG1.hic.250624.v3mitohifi"])
            entry = grouped["OG1.hic.250624.v3mitohifi"]
            self.assertEqual(
                set(entry["steps"]), {"assembly", "ena_validation", "qc_validator"}
            )
            self.assertEqual(entry["annotation"], "emma102")

            output = root / "summary.tsv"
            MODULE.write_summary_tsv(grouped, output)
            header, row = (line.split("\t") for line in output.read_text().splitlines())
            values = dict(zip(header, row))
            self.assertEqual(values["ena_validation"], "success")
            self.assertEqual(values["qc_validator"], "validator_2_set")


if __name__ == "__main__":
    unittest.main()
