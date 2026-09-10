"""Tests for the lca_validation writer.

This behaviour used to live inside species_validation.py, where it was welded to
the species comparison that decides the QC gate. Splitting it out is what lets that
comparison run with no database; these tests pin the parts that moved, so the split
cannot quietly drop the overwrite guard along the way.
"""

import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "push_species_validation", ROOT / "bin" / "push_species_validation.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


KEY = ("OG470", "ilmn", "230607", "getorg1770", "emma102")


class FakeCursor:
    def __init__(self, existing_row):
        self.existing_row = existing_row
        self.queries = []
        self.params = []

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, params=None):
        self.queries.append(query)
        self.params.append(params)

    def fetchone(self):
        return self.existing_row


class FakeConnection:
    def __init__(self, existing_row):
        self.cursor_obj = FakeCursor(existing_row)
        self.committed = False
        self.rolled_back = False

    def cursor(self):
        return self.cursor_obj

    def commit(self):
        self.committed = True

    def rollback(self):
        self.rolled_back = True

    def close(self):
        pass


class UpsertTests(unittest.TestCase):
    def _run(self, existing_row=None, force=False, **kwargs):
        conn = FakeConnection(existing_row)

        class FakePsycopg2:
            connect = staticmethod(lambda **_kw: conn)

        original = MODULE.psycopg2
        MODULE.psycopg2 = FakePsycopg2
        try:
            buffer = io.StringIO()
            with redirect_stdout(buffer):
                ok = MODULE.upsert_lca_validation(
                    {"dbname": "t"}, KEY, kwargs.pop("species", "Genus species"),
                    force=force, **kwargs
                )
        finally:
            MODULE.psycopg2 = original
        return ok, conn, buffer.getvalue()

    def test_it_writes_the_row_and_commits(self):
        ok, conn, _out = self._run(validated_rank="species")
        self.assertTrue(ok)
        self.assertTrue(conn.committed)
        inserts = [q for q in conn.cursor_obj.queries if "INSERT INTO lca_validation" in q]
        self.assertEqual(len(inserts), 1)
        params = conn.cursor_obj.params[-1]
        self.assertEqual(params["og_id"], "OG470")
        self.assertEqual(params["annotation"], "emma102")
        self.assertEqual(params["validated_species_name"], "Genus species")
        self.assertEqual(params["validated_rank"], "species")

    def test_a_human_validator_is_preserved(self):
        """The guard that matters: an automated re-run must not overwrite a row a
        reviewer set by hand."""
        ok, conn, out = self._run(existing_row=("curator", "Real species"))
        self.assertTrue(ok)
        self.assertFalse(conn.committed)
        self.assertFalse(any("INSERT INTO" in q for q in conn.cursor_obj.queries))
        self.assertIn("Existing values preserved", out)

    def test_force_overwrites_a_human_validator(self):
        ok, conn, out = self._run(existing_row=("curator", "Real species"), force=True)
        self.assertTrue(ok)
        self.assertTrue(conn.committed)
        self.assertTrue(any("INSERT INTO" in q for q in conn.cursor_obj.queries))
        self.assertIn("--force", out)

    def test_an_existing_nf_core_row_is_overwritten_without_force(self):
        ok, conn, _out = self._run(existing_row=("nf-core", "Old species"))
        self.assertTrue(ok)
        self.assertTrue(conn.committed)
        self.assertTrue(any("INSERT INTO" in q for q in conn.cursor_obj.queries))

    def test_a_failed_write_returns_false(self):
        """Must propagate into a non-zero exit: a failed write that returned normally
        is how rows lost to the FK write-ordering race stayed invisible."""
        class ExplodingPsycopg2:
            @staticmethod
            def connect(**_kw):
                raise RuntimeError("connection refused")

        original = MODULE.psycopg2
        MODULE.psycopg2 = ExplodingPsycopg2
        try:
            with redirect_stdout(io.StringIO()):
                ok = MODULE.upsert_lca_validation({"dbname": "t"}, KEY, "Genus species")
        finally:
            MODULE.psycopg2 = original
        self.assertFalse(ok)


class RecordDispatchTests(unittest.TestCase):
    """main() reads the record species_validation.py wrote and does what it says."""

    def _main(self, record, force=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            record_path = root / "validation_record.json"
            record_path.write_text(json.dumps(record))
            config_path = root / "db.cfg"
            config_path.write_text(
                "[postgres]\ndbname=t\nuser=u\npassword=p\nhost=h\nport=5432\n"
            )

            calls = []

            def fake_upsert(*args, **kwargs):
                calls.append((args, kwargs))
                return True

            original_upsert = MODULE.upsert_lca_validation
            original_argv = sys.argv
            MODULE.upsert_lca_validation = fake_upsert
            sys.argv = ["push_species_validation.py"] + (["--force"] if force else []) + [
                str(config_path), str(record_path)
            ]
            try:
                buffer = io.StringIO()
                with redirect_stdout(buffer):
                    rc = MODULE.main()
            finally:
                MODULE.upsert_lca_validation = original_upsert
                sys.argv = original_argv
            return rc, calls, buffer.getvalue()

    def test_a_skip_record_writes_nothing_and_succeeds(self):
        """A sample the gate held is not an error, and must leave the table alone."""
        rc, calls, out = self._main(
            {"action": "skip", "reason": "not validated (no Found_in_blast_YN=Yes)"}
        )
        self.assertEqual(rc, 0)
        self.assertEqual(calls, [])
        self.assertIn("No lca_validation row to write", out)

    def test_an_upsert_record_is_passed_through_field_for_field(self):
        rc, calls, _out = self._main({
            "action": "upsert",
            "key": {"og_id": "OG470", "tech": "ilmn", "seq_date": "230607",
                    "code": "getorg1770", "annotation": "emma102"},
            "validated_species_name": "Genus species",
            "validator": "nf-core",
            "validated_rank": "family",
            "lca_genus": "Lamprogrammus",
        })
        self.assertEqual(rc, 0)
        self.assertEqual(len(calls), 1)
        args, kwargs = calls[0]
        self.assertEqual(args[1], KEY)
        self.assertEqual(args[2], "Genus species")
        self.assertEqual(kwargs["validated_rank"], "family")
        self.assertEqual(kwargs["lca_genus"], "Lamprogrammus")

    def test_force_reaches_the_upsert(self):
        _rc, calls, _out = self._main({
            "action": "upsert",
            "key": {"og_id": "OG470", "tech": "ilmn", "seq_date": "230607",
                    "code": "getorg1770", "annotation": "emma102"},
            "validated_species_name": None,
        }, force=True)
        self.assertTrue(calls[0][1]["force"])

    def test_an_upsert_with_no_key_fails_loudly(self):
        rc, calls, _out = self._main({"action": "upsert", "key": None})
        self.assertEqual(rc, 1)
        self.assertEqual(calls, [])

    def test_an_unknown_action_fails_rather_than_guessing(self):
        rc, calls, _out = self._main({"action": "maybe"})
        self.assertEqual(rc, 1)
        self.assertEqual(calls, [])


if __name__ == "__main__":
    unittest.main()
