"""A database upload that fails must fail its task.

push_mtdna_assm_results.py and push_emma_annotation_results.py used to print
"❌ Database error" and exit 0, so Nextflow recorded a failed upload as a success,
cached it, and -resume never re-ran it (genomes_db_review docs/r6_merge_record.md: 20
uploads of RELEASE_TEST_261102 lost that way). They now exit 1. The run still carries
on: the default errorStrategy retries the task once and then ignores it, and the
completion handler lists it in sql_uploaded_data/FAILED_UPLOADS.tsv.

The database error here is a connection refused on a closed port, so these run the real
scripts end to end. They need pandas, numpy and psycopg2, and skip without them.
"""

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "unit" / "fixtures" / "upload_exit"
HAVE_DEPS = all(importlib.util.find_spec(m) for m in ("pandas", "numpy", "psycopg2"))

DEAD_CFG = "[postgres]\ndbname = x\nuser = x\npassword = x\nhost = 127.0.0.1\nport = 1\n"


@unittest.skipUnless(HAVE_DEPS, "needs pandas, numpy and psycopg2")
class DatabaseErrorExitsNonZero(unittest.TestCase):
    def run_script(self, script, *args):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "dead.cfg").write_text(DEAD_CFG)
            (tmp / "a.fasta").write_text(">OG1.ilmn.250101.getorg1770\nACGTACGTACGT\n")
            (tmp / "get_org.log.txt").write_text("Result status of animal_mt: circular genome\n")
            result = subprocess.run(
                [sys.executable, str(ROOT / "bin" / script), *[a.format(tmp=tmp) for a in args]],
                cwd=tmp, capture_output=True, text=True, timeout=120)
        return result

    def test_assembly_upload(self):
        result = self.run_script("push_mtdna_assm_results.py", "{tmp}/dead.cfg",
                                 "OG1.ilmn.250101.getorg1770", "{tmp}/get_org.log.txt", "{tmp}/a.fasta")
        self.assertIn("Database error", result.stdout)
        self.assertEqual(result.returncode, 1, result.stdout[-500:])

    def test_annotation_upload(self):
        stats = FIXTURES / "OG3150.ilmn.260925.getorg1770.annotation_stats.csv"
        result = self.run_script("push_emma_annotation_results.py", "{tmp}/dead.cfg", "OG3150", str(stats))
        self.assertIn("Database error", result.stdout)
        self.assertEqual(result.returncode, 1, result.stdout[-500:])


class ModulesPassTheExitCodeOn(unittest.TestCase):
    """The two modules that forced `exit 0` now exit with the script's code."""

    def test_no_forced_success(self):
        for module in ("qc_validator", "ena_validation"):
            text = (ROOT / "modules" / "local" / "upload_results" / module / "main.nf").read_text()
            script = text.split("script:", 1)[1].split("stub:", 1)[0]
            self.assertNotIn("\n    exit 0\n", script, module)
            self.assertIn("exit \\${upload_rc}", script, module)


if __name__ == "__main__":
    unittest.main()
