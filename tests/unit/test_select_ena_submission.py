import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "select_ena_submission", ROOT / "bin" / "select_ena_submission.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def package(seqid, digest, status="READY", local="PASS"):
    return {
        "og_id": seqid.split(".")[0],
        "full_seqid": seqid,
        "normalised_circular_sha256": digest,
        "package_status": status,
        "local_validation_status": local,
    }


def write_metadata(directory, seqid, extra=None):
    directory.mkdir(parents=True, exist_ok=True)
    payload = package(seqid, "a" * 64)
    payload.update(extra or {})
    path = directory / f"{seqid}.package_metadata.json"
    path.write_text(json.dumps(payload))
    return path


class PackagePathTests(unittest.TestCase):
    """Where the database is told the package lives.

    Nextflow stages metadata JSONs flat into a task work directory, so the
    staged parent is never the answer: it is thrown away at the end of the run,
    and the ENA submitter reads package_path out of the database days or weeks
    later.
    """

    def test_published_path_beats_the_staged_location(self):
        seqid = "OG910.hifi.250101.v3mitohifi.emma102"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            published = root / "results" / "ena" / "package"
            published.mkdir(parents=True)
            staged = write_metadata(
                root / "work" / "ab" / "cdef",
                seqid,
                {"published_package_path": str(published)},
            )
            # Webin results sit beside the published package, never beside the
            # staged copy. Reading them from the staged parent is what let a
            # NOT_RUN overwrite a recorded PASS.
            validation = root / "results" / "ena" / "validation" / "webin_test"
            validation.mkdir(parents=True)
            (validation / f"{seqid}.webin_test_status.tsv").write_text(
                f"full_seqid\tstatus\n{seqid}\tPASS\n"
            )
            loaded = MODULE.load_packages([staged])[0]
            self.assertEqual(loaded["_package_path"], str(published))
            self.assertEqual(loaded["_metadata_path"], str(staged))
            self.assertEqual(loaded["webin_test_status"], "PASS")

    def test_missing_published_path_falls_back_and_warns(self):
        seqid = "OG5.ilmn.260101.getorg1770.emma102"
        with tempfile.TemporaryDirectory() as tmp:
            staged = write_metadata(Path(tmp) / "package", seqid)
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                loaded = MODULE.load_packages([staged])[0]
            self.assertEqual(loaded["_package_path"], str(staged.parent))
            self.assertIn("no published_package_path", stderr.getvalue())

    def test_absent_webin_results_leave_the_status_unset(self):
        # load_packages must not invent NOT_RUN; the caller has to be able to
        # tell "no news" from "reset to NOT_RUN".
        seqid = "OG5.ilmn.260101.getorg1770.emma102"
        with tempfile.TemporaryDirectory() as tmp:
            staged = write_metadata(Path(tmp) / "package", seqid)
            with contextlib.redirect_stderr(io.StringIO()):
                loaded = MODULE.load_packages([staged])[0]
            self.assertNotIn("webin_test_status", loaded)
            self.assertNotIn("webin_production_status", loaded)


class SelectEnaSubmissionTests(unittest.TestCase):
    def test_every_technology_is_selected_independently(self):
        packages = [
            package("OG910.ilmn.260101.getorg1770.emma102", "a" * 64),
            package("OG910.hifi.250101.v3mitohifi.emma102", "b" * 64),
            package("OG910.hic.250101.getorg1770.emma102", "c" * 64),
        ]
        _, results = MODULE.build_report(packages, {})
        self.assertEqual(
            sorted(results), [("OG910", "hic"), ("OG910", "hifi"), ("OG910", "ilmn")]
        )
        for key, result in results.items():
            self.assertEqual(result["status"], "SELECTED", key)
        self.assertEqual(
            results[("OG910", "hifi")]["selected"],
            "OG910.hifi.250101.v3mitohifi.emma102",
        )

    def test_identical_sequence_in_two_technologies_is_published_twice(self):
        # Both records are genuinely separate assemblies from separate data, so
        # neither suppresses the other even when the sequences are identical.
        packages = [
            package("OG910.hifi.250101.v3mitohifi.emma102", "a" * 64),
            package("OG910.hic.250101.getorg1770.emma102", "a" * 64),
        ]
        _, results = MODULE.build_report(packages, {})
        self.assertEqual(len(results), 2)
        self.assertEqual(
            {result["selected"] for result in results.values()},
            {seqid["full_seqid"] for seqid in packages},
        )

    def test_one_selection_per_specimen_and_technology(self):
        packages = [
            package("OG910.hifi.250101.v3mitohifi.emma102", "a" * 64),
            package("OG910.hifi.250505.v323mitohifi.emma102", "a" * 64),
            package("OG910.hifi.240101.v321mitohifi.emma102", "a" * 64),
        ]
        report, results = MODULE.build_report(packages, {})
        self.assertEqual(list(results), [("OG910", "hifi")])
        # Equivalent sequences: newest sequencing date wins the tiebreak.
        self.assertEqual(
            results[("OG910", "hifi")]["selected"],
            "OG910.hifi.250505.v323mitohifi.emma102",
        )
        self.assertEqual(
            results[("OG910", "hifi")]["reason"],
            "equivalent_candidates_date_seqid_tiebreak",
        )
        self.assertEqual(sum(row["selected"] == "true" for row in report), 1)

    def test_distinct_sequences_within_one_technology_require_review(self):
        result = MODULE.select_group(
            [
                package("OG1.hifi.260101.v3mitohifi.emma102", "a" * 64),
                package("OG1.hifi.250101.v321mitohifi.emma102", "b" * 64),
            ]
        )
        self.assertEqual(result["status"], "MANUAL_REVIEW_REQUIRED")
        self.assertIsNone(result["selected"])

    def test_only_ready_local_pass_candidates_are_eligible(self):
        result = MODULE.select_group(
            [
                package("OG1.hifi.260101.v3mitohifi.emma102", "a" * 64, local="FAIL"),
                package("OG1.hifi.250101.v321mitohifi.emma102", "b" * 64),
            ]
        )
        self.assertEqual(result["selected"], "OG1.hifi.250101.v321mitohifi.emma102")

    def test_unknown_technology_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "expected one of"):
            MODULE.tech_of("OG1.ont.260101.flye.emma102")

    def test_manual_decision_applies_to_one_technology_only(self):
        packages = [
            package("OG1.hifi.260101.v3mitohifi.emma102", "a" * 64),
            package("OG1.hifi.250101.v321mitohifi.emma102", "b" * 64),
            package("OG1.ilmn.260101.getorg1770.emma102", "c" * 64),
        ]
        report, results = MODULE.build_report(
            packages,
            {
                ("OG1", "hifi"): {
                    "og_id": "OG1",
                    "selected_seqid": packages[1]["full_seqid"],
                    "reviewer": "curator",
                    "reason": "read support",
                }
            },
        )
        self.assertEqual(results[("OG1", "hifi")]["selected"], packages[1]["full_seqid"])
        self.assertEqual(results[("OG1", "hifi")]["reason"], "manual:read support")
        # The Illumina group is untouched by the HiFi decision.
        self.assertEqual(
            results[("OG1", "ilmn")]["selected"], packages[2]["full_seqid"]
        )
        self.assertEqual(sum(row["selected"] == "true" for row in report), 2)

    def test_decision_file_rejects_two_rows_for_one_technology(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "decisions.tsv"
            path.write_text(
                "og_id\tselected_seqid\treviewer\treason\n"
                "OG1\tOG1.hifi.260101.v3mitohifi.emma102\tcurator\tfirst\n"
                "OG1\tOG1.hifi.250101.v321mitohifi.emma102\tcurator\tsecond\n"
            )
            with self.assertRaisesRegex(ValueError, "Duplicate manual decision"):
                MODULE.read_decisions(path)

    def test_decision_file_accepts_one_row_per_technology(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "decisions.tsv"
            path.write_text(
                "og_id\tselected_seqid\treviewer\treason\n"
                "OG1\tOG1.hifi.260101.v3mitohifi.emma102\tcurator\thifi pick\n"
                "OG1\tOG1.ilmn.260101.getorg1770.emma102\tcurator\tilmn pick\n"
            )
            decisions = MODULE.read_decisions(path)
            self.assertEqual(sorted(decisions), [("OG1", "hifi"), ("OG1", "ilmn")])

    def test_report_digest_is_stable(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.tsv"
            rows, _ = MODULE.build_report(
                [package("OG1.hifi.260101.v3mitohifi.emma102", "a" * 64)], {}
            )
            first = MODULE.write_report(path, rows)
            second = MODULE.write_report(path, rows)
            self.assertEqual(first, second)

    def test_selected_package_sheet_points_to_finished_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            package_dir = Path(tmp) / "package"
            package_dir.mkdir()
            metadata = package_dir / "candidate.package_metadata.json"
            candidate = package("OG910.hifi.250101.v3mitohifi.emma102", "a" * 64)
            candidate.update(
                {
                    "_package_path": str(package_dir),
                    "_metadata_path": str(metadata),
                    "study": "PRJEB123419",
                    "biosample_accession": "SAMEA1",
                }
            )
            _, results = MODULE.build_report([candidate], {})
            selected = Path(tmp) / "selected.tsv"
            MODULE.write_selected_packages(selected, [candidate], results)
            with selected.open() as handle:
                row = list(MODULE.csv.DictReader(handle, delimiter="\t"))[0]
            self.assertEqual(row["package_path"], str(package_dir))
            self.assertEqual(row["full_seqid"], candidate["full_seqid"])
            self.assertEqual(row["tech"], "hifi")
            self.assertEqual(row["study"], "PRJEB123419")


if __name__ == "__main__":
    unittest.main()
