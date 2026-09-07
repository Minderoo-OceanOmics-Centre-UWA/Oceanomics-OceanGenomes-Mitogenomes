"""One unrepresentable row must not cost a sample its whole upload.

The regression these cover lost 87 assemblies their LCA rows across
batch-12 .. batch-20: the push scripts caught per-row exceptions, but without a
savepoint PostgreSQL had already aborted the transaction, so every later row
failed with "current transaction is aborted" and the commit was downgraded to a
rollback -- while the script printed a tally of "succeeded" rows and exited 0.

FakeConnection reproduces exactly that behaviour, so a script without savepoints
fails these tests the way the database did.
"""

import importlib.util
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Unlike most of this suite these two scripts are NOT stdlib-only: both read their
# input with pandas and both import psycopg2 at module scope, and the functions
# under test call pd.read_csv directly, so stubbing the imports away would not let
# them run. Guard the load instead, the same way test_mitos_to_emma.py guards
# Biopython. Without this the module raises at import and the whole file is
# reported as an ERROR -- which reads as a broken suite rather than as a missing
# dependency, and masks any real failure alongside it.
try:
    import pandas  # noqa: F401
    import psycopg2  # noqa: F401
    DEPS_AVAILABLE = True
except ImportError:
    DEPS_AVAILABLE = False

if DEPS_AVAILABLE:
    BLAST = load("push_lca_blast_results")
    RAW = load("push_lca_raw_results")


class AbortedTransaction(Exception):
    pass


class FakeCursor:
    """A cursor with PostgreSQL's actual failure semantics.

    Once a statement raises, the transaction is aborted and every subsequent
    statement raises too -- until a ROLLBACK TO SAVEPOINT clears it.
    """

    def __init__(self, connection, reject):
        self.connection = connection
        self.reject = reject

    def execute(self, query, params=None):
        stripped = query.strip()
        if stripped.startswith("ROLLBACK TO SAVEPOINT"):
            self.connection.aborted = False
            # Everything written since the savepoint is discarded.
            self.connection.pending.pop()
            return
        if self.connection.aborted:
            raise AbortedTransaction(
                "current transaction is aborted, commands ignored until end of "
                "transaction block"
            )
        if stripped.startswith("SAVEPOINT"):
            self.connection.pending.append([])
            return
        if stripped.startswith("RELEASE SAVEPOINT"):
            released = self.connection.pending.pop()
            self.connection.pending[-1].extend(released)
            return
        if self.reject(params):
            self.connection.aborted = True
            raise AbortedTransaction('invalid input syntax for type integer: "1;2"')
        self.connection.pending[-1].append(params)

    def fetchone(self):
        return (True,)

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False


class FakeConnection:
    """Commits the pending statements, unless the transaction was aborted."""

    def __init__(self, reject):
        self.reject = reject
        self.aborted = False
        # A stack of statement lists: the outermost transaction, then one per
        # open savepoint.
        self.pending = [[]]
        self.committed = []

    def cursor(self):
        return FakeCursor(self, self.reject)

    def commit(self):
        if self.aborted:
            # What PostgreSQL does: COMMIT on an aborted transaction rolls back.
            self.pending = [[]]
            return
        self.committed.extend(self.pending[0])
        self.pending = [[]]

    def rollback(self):
        self.pending = [[]]

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *_args):
        if exc_type is None:
            self.commit()
        else:
            self.rollback()
        return False


BLAST_ROWS = [
    # query_id, match_sequence_id, taxon_id, then 21 more columns.
    ("OG1.ilmn.260514.getorg1770.mitos2110", "gb|A1|", "111"),
    ("OG1.ilmn.260514.getorg1770.mitos2110", "gb|A2|", "222;333"),  # the poison row
    ("OG1.ilmn.260514.getorg1770.mitos2110", "gb|A3|", "444"),
]


def write_blast_tsv(path):
    with open(path, "w") as handle:
        for query_id, accession, taxon_id in BLAST_ROWS:
            filler = ["1"] * (len(BLAST.blast_column_headers) - 5)
            handle.write(
                "\t".join([query_id, accession, taxon_id, "Genus species", "common"]
                          + filler[:-1] + ["12s"])
                + "\n"
            )


@unittest.skipUnless(DEPS_AVAILABLE, "pandas/psycopg2 not installed")
class BlastRowIsolationTests(unittest.TestCase):
    def run_push(self, reject):
        connection = FakeConnection(reject)
        original = BLAST.psycopg2.connect
        BLAST.psycopg2.connect = lambda **_kwargs: connection
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tsv = Path(tmp) / "blast_combined.tsv"
                write_blast_tsv(tsv)
                out = io.StringIO()
                with redirect_stdout(out):
                    failures = BLAST.process_blast(str(tsv), "OG1", {})
            return connection, failures, out.getvalue()
        finally:
            BLAST.psycopg2.connect = original

    def test_good_rows_commit_despite_a_rejected_row(self):
        connection, failures, output = self.run_push(
            lambda params: ";" in str(params["taxon_id"])
        )
        self.assertEqual(failures, 1)
        self.assertEqual(len(connection.committed), 2)
        self.assertNotIn("gb|A2|", [p["match_sequence_id"] for p in connection.committed])
        self.assertIn("finished with errors", output)
        # The tick is reserved for an upload that actually landed in full.
        self.assertNotIn("✅", output)

    def test_clean_upload_commits_everything(self):
        connection, failures, output = self.run_push(lambda _params: False)
        self.assertEqual(failures, 0)
        self.assertEqual(len(connection.committed), 3)
        self.assertIn("✅", output)


RAW_HEADER = [
    "seq_id", "accession_id", "sequence_region", "lca_run_date", "confidence_score",
]


def write_raw_tsv(path, accessions):
    with open(path, "w") as handle:
        handle.write("\t".join(RAW_HEADER) + "\n")
        for accession in accessions:
            handle.write(
                "\t".join(
                    ["OG1.ilmn.260514.getorg1770.mitos2110", accession, "12s",
                     "260902", "0.0"]
                )
                + "\n"
            )


@unittest.skipUnless(DEPS_AVAILABLE, "pandas/psycopg2 not installed")
class RawRowIsolationTests(unittest.TestCase):
    def test_good_rows_commit_despite_a_rejected_row(self):
        connection = FakeConnection(lambda params: params["accession_id"] == "gb|B2|")
        with tempfile.TemporaryDirectory() as tmp:
            tsv = Path(tmp) / "lca_raw.12s.tsv"
            write_raw_tsv(tsv, ["gb|B1|", "gb|B2|", "gb|B3|"])
            written = {}
            out = io.StringIO()
            with redirect_stdout(out):
                cursor = connection.cursor()
                inserted, refreshed, failed = RAW.process_file(cursor, str(tsv), written)
            connection.commit()

        self.assertEqual((inserted, refreshed, failed), (2, 0, 1))
        self.assertEqual(
            [p["accession_id"] for p in connection.committed], ["gb|B1|", "gb|B3|"]
        )


if __name__ == "__main__":
    unittest.main()
