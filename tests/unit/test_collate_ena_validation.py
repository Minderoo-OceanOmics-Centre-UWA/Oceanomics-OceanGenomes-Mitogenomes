import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "collate_ena_validation", ROOT / "bin" / "collate_ena_validation.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_tsv(path, row):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    return path


class CollateEnaValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.prefix = "OG123.hifi.260101.final"

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, paths, **overrides):
        settings = dict(
            assembly_prefix=self.prefix,
            og_id="OG123",
            ena_study="PRJEB1",
            validation_mode="pipeline",
            validation_attempt="attempt1",
            webin_requested=True,
        )
        settings.update(overrides)
        return MODULE.build_record(paths, **settings)

    def table(self, status="PASS", **counts):
        row = dict(
            sample=self.prefix, status=status, reject_count="0", error_count="0",
            warning_count="0", info_count="0", fatal_discrepancy_count="0",
            nostop_count="0", blocking_codes="", warning_codes="",
        )
        row.update({key: str(value) for key, value in counts.items()})
        return write_tsv(self.root / f"{self.prefix}.table2asn_status.tsv", row)

    def conversion(self, status="PASS", reason="ok", exit_code="0"):
        return write_tsv(
            self.root / f"{self.prefix}.ena_conversion_status.tsv",
            dict(sample=self.prefix, status=status, reason=reason, seqret_exit=exit_code),
        )

    def webin(self, status="PASS", reason="validated", exit_code="0"):
        return write_tsv(
            self.root / f"{self.prefix}.webin_status.tsv",
            dict(sample=self.prefix, status=status, reason=reason, webin_exit=exit_code),
        )

    def test_pass_record_and_artifact_provenance(self):
        flatfile = self.root / f"{self.prefix}.embl.gz"
        flatfile.write_bytes(b"flatfile")
        manifest = self.root / f"{self.prefix}.webin_manifest.txt"
        manifest.write_text("STUDY\tPRJEB1\n")
        record = self.build([self.table(), self.conversion(), self.webin(), flatfile, manifest])
        self.assertEqual(record["submission_ready"], "true")
        self.assertEqual(record["webin_status"], "PASS")
        self.assertEqual(len(record["flatfile_sha256"]), 64)
        self.assertNotIn(str(self.root), record["flatfile_name"])

    def test_table2asn_failure_skips_downstream(self):
        record = self.build([self.table("FAIL_TABLE2ASN", error_count=1, blocking_codes="SEQ_FEAT.NoStop")])
        self.assertEqual(record["conversion_status"], "SKIPPED_TABLE2ASN")
        self.assertEqual(record["webin_status"], "NOT_RUN")
        self.assertEqual(record["nostop_count"], "0")
        self.assertEqual(record["submission_ready"], "false")

    def test_conversion_failure(self):
        record = self.build([self.table(), self.conversion("FAIL_CONVERSION", "seqret_exit_2", "2")])
        self.assertEqual(record["conversion_status"], "FAIL_CONVERSION")
        self.assertEqual(record["conversion_exit"], "2")
        self.assertEqual(record["webin_status"], "NOT_RUN")

    def test_warning_only_table2asn_continues(self):
        record = self.build([
            self.table("PASS", warning_count=2, warning_codes="SEQ_DESCR.SuspiciousComplete"),
            self.conversion(),
        ])
        self.assertEqual(record["table2asn_status"], "PASS")
        self.assertEqual(record["warning_count"], "2")
        self.assertEqual(record["conversion_status"], "PASS")

    def test_webin_validation_and_infrastructure_failures(self):
        for status, reason in (("FAIL_WEBIN", "validation_failed"), ("FAIL_INFRASTRUCTURE", "network")):
            with self.subTest(status=status):
                record = self.build([self.table(), self.conversion(), self.webin(status, reason, "1")])
                self.assertEqual(record["webin_status"], status)
                self.assertEqual(record["webin_reason"], reason)
                self.assertEqual(record["submission_ready"], "false")

    def test_webin_disabled_is_explicit(self):
        record = self.build([self.table(), self.conversion()], webin_requested=False)
        self.assertEqual(record["webin_status"], "NOT_REQUESTED")
        self.assertEqual(record["webin_reason"], "not_requested")

    def test_validate_mode_preflight_failure(self):
        status = write_tsv(
            self.root / f"{self.prefix}.ena_preflight_status.tsv",
            dict(sample=self.prefix, status="FAIL_PREFLIGHT", reason="invalid_gzip"),
        )
        checks = write_tsv(
            self.root / f"{self.prefix}.ena_preflight_check.tsv",
            dict(sample=self.prefix, gzip_exit="1"),
        )
        record = self.build([status, checks], validation_mode="validate")
        self.assertEqual(record["table2asn_status"], "NOT_APPLICABLE")
        self.assertEqual(record["conversion_status"], "NOT_APPLICABLE")
        self.assertEqual(record["preflight_status"], "FAIL_PREFLIGHT")
        self.assertEqual(record["preflight_exit"], "1")

    def test_malformed_table_status_is_retained(self):
        malformed = write_tsv(
            self.root / f"{self.prefix}.table2asn_status.tsv", dict(sample=self.prefix, reason="bad")
        )
        record = self.build([malformed])
        self.assertEqual(record["table2asn_status"], "MALFORMED_STATUS")
        self.assertEqual(record["conversion_status"], "SKIPPED_TABLE2ASN")

    def test_digest_ignores_workflow_provenance_but_changes_with_result(self):
        paths = [self.table(), self.conversion(), self.webin()]
        first = self.build(paths, workflow_run_name="run1")
        second = self.build(paths, workflow_run_name="run2", workflow_session_id="new")
        failed = self.build([paths[0], paths[1], self.webin("FAIL_WEBIN", "bad", "1")])
        self.assertEqual(first["result_digest"], second["result_digest"])
        self.assertNotEqual(first["result_digest"], failed["result_digest"])

    def test_batch_summary_keeps_multiple_assemblies(self):
        rows = []
        for suffix in ("a", "b"):
            record = self.build([self.table()], assembly_prefix=f"OG123.hifi.260101.{suffix}")
            path = self.root / f"{suffix}.ena_validation_result.tsv"
            MODULE.write_rows(path, MODULE.RECORD_COLUMNS, [record])
            rows.append(path)
        read = MODULE.read_records(rows)
        self.assertEqual([row["assembly_prefix"] for row in read], ["OG123.hifi.260101.a", "OG123.hifi.260101.b"])
        self.assertEqual(MODULE.expand_patterns([str(self.root / "*.ena_validation_result.tsv")]), rows)


if __name__ == "__main__":
    unittest.main()
