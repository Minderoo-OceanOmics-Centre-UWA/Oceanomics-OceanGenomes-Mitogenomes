"""The overwrite policy for an existing mitogenome_data row.

`params.force_db_overwrite` defaults to true, so these cases are what stands
between a re-run and the production record. The two unconditional refusals below
are the point of the test: they must hold with force ON, because a re-run whose
assembly happened to fail must never replace a good result with 'failed to
assemble'.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))


def _load_with_stubbed_deps(name):
    """Import a bin/ script that needs pandas/psycopg2, without either installed.

    decide_write() is stdlib-only; the third-party imports it is sitting behind are
    module-scope and never reached by the function under test, so empty stubs are
    enough to get the module loaded. This is the opposite call from
    test_push_lca_row_isolation.py, which skips instead because the functions IT
    tests really do call pd.read_csv.

    sys.modules is restored afterwards so a stub cannot leak into another test
    module and convince it the real library is present.
    """
    stubbed = [dep for dep in ("pandas", "numpy", "psycopg2") if dep not in sys.modules]
    for dep in stubbed:
        sys.modules[dep] = types.ModuleType(dep)
    try:
        spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        for dep in stubbed:
            del sys.modules[dep]


MODULE = _load_with_stubbed_deps("push_mtdna_assm_results")

decide_write = MODULE.decide_write
FAILED = "failed to assemble"
GOOD = "circular genome"


class DecideWriteTests(unittest.TestCase):
    def test_successful_row_is_never_downgraded_to_failure(self):
        """The guard that force must not be able to switch off."""
        for force in (False, True):
            with self.subTest(force=force):
                should_write, reason = decide_write(GOOD, True, force)
                self.assertFalse(should_write)
                self.assertEqual(reason, MODULE.SKIP_WOULD_DOWNGRADE)

    def test_failed_row_plus_failed_run_is_a_no_op_either_way(self):
        for force in (False, True):
            with self.subTest(force=force):
                should_write, reason = decide_write(FAILED, True, force)
                self.assertFalse(should_write)
                self.assertEqual(reason, MODULE.SKIP_ALREADY_FAILED)

    def test_case_and_whitespace_variants_still_read_as_failed(self):
        for stored in ("Failed To Assemble", "  failed to assemble  ", "FAILED TO ASSEMBLE"):
            with self.subTest(stored=stored):
                self.assertEqual(decide_write(stored, True, True)[1], MODULE.SKIP_ALREADY_FAILED)
                # ... and a success still upgrades it.
                self.assertTrue(decide_write(stored, False, False)[0])

    def test_success_upgrades_a_failed_row_without_force(self):
        should_write, reason = decide_write(FAILED, False, False)
        self.assertTrue(should_write)
        self.assertIsNone(reason)

    def test_incomplete_row_is_always_completed(self):
        for failed_to_assemble in (False, True):
            for force in (False, True):
                with self.subTest(failed=failed_to_assemble, force=force):
                    should_write, reason = decide_write(None, failed_to_assemble, force)
                    self.assertTrue(should_write)
                    self.assertIsNone(reason)

    def test_force_is_what_unlocks_overwriting_a_real_result(self):
        should_write, reason = decide_write(GOOD, False, False)
        self.assertFalse(should_write)
        self.assertEqual(reason, MODULE.SKIP_ROW_EXISTS)

        should_write, reason = decide_write(GOOD, False, True)
        self.assertTrue(should_write)
        self.assertIsNone(reason)


class SkipReasonContractTests(unittest.TestCase):
    """compile_upload_report.py classifies these rows by substring match."""

    def setUp(self):
        spec = importlib.util.spec_from_file_location(
            "compile_upload_report", ROOT / "bin" / "compile_upload_report.py"
        )
        self.report = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.report)

    def _classify(self, reason):
        body = f"⚠️ Existing values preserved for OG123.ilmn.240101.0: {reason}."
        return self.report.classify_status("assembly", body)

    def test_downgrade_refusal_is_reported_as_a_kept_prior_success(self):
        self.assertEqual(
            self._classify(MODULE.SKIP_WOULD_DOWNGRADE), "failed_kept_prior_success"
        )

    def test_row_exists_skip_is_still_recognised(self):
        self.assertNotEqual(self._classify(MODULE.SKIP_ROW_EXISTS), "unknown")


if __name__ == "__main__":
    unittest.main()
