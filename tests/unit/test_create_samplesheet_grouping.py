"""Each library must become its own samplesheet row, not be merged by sample.

group_files_by_sample() used to key only on the sample id, so every library a
specimen had collapsed into a single row and CAT_FASTQ concatenated them into
one assembly. Two archived cases hit this: specimens with two ilmn dates, and
Hi-C tubes sequenced on more than one run, where query_hic_date() can only
return the most recent seq_date for the tube regardless of which run a file
came from.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"


def load_create_samplesheet():
    """Import create_samplesheet.py without requiring psycopg2 to be installed."""
    spec = importlib.util.spec_from_file_location(
        "create_samplesheet_under_test", BIN / "create_samplesheet.py"
    )
    module = importlib.util.module_from_spec(spec)
    patches = {}
    try:
        import psycopg2  # noqa: F401
    except ImportError:
        patches["psycopg2"] = types.ModuleType("psycopg2")
    sys.path.insert(0, str(BIN))          # for species_name_utils
    try:
        with mock.patch.dict(sys.modules, patches):
            spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(BIN))
    return module


cs = load_create_samplesheet()


class HicRunDateTests(unittest.TestCase):
    def test_run_date_after_hicl_token_is_found(self):
        self.assertEqual(
            cs.extract_hic_run_date("OG2088G-1_HICL_260605_S3_L001_R1_001.fastq.gz"),
            "260605",
        )

    def test_absent_run_date_returns_empty(self):
        # Freshly sequenced reads keep sequencer naming and fall back to the DB.
        self.assertEqual(
            cs.extract_hic_run_date("OG2088G-1_HICL_S3_L001_R1_001.fastq.gz"), ""
        )

    def test_lane_and_sample_fields_are_not_mistaken_for_a_date(self):
        self.assertEqual(
            cs.extract_hic_run_date("OG820M-1_HICL_S4_L001_R1_001.fastq.gz"), ""
        )

    def test_sample_id_is_unchanged_by_the_inserted_date(self):
        # The DB tube lookup keys on this, so it must stay OG<n><tube>_HICL.
        self.assertEqual(
            cs.extract_sample_info("OG2088G-1_HICL_260605_S3_L001_R1_001.fastq.gz"),
            "OG2088G-1_HICL",
        )

    def test_technology_still_parses_as_hic(self):
        self.assertEqual(
            cs.determine_sequencing_type("OG2088G-1_HICL_260605_S3_L001_R1_001.fastq.gz"),
            "hic",
        )

    def test_read_direction_still_parses(self):
        self.assertEqual(
            cs.determine_file_type("OG2088G-1_HICL_260605_S3_L001_R2_001.fastq.gz"), "R2"
        )


class LibraryGroupingTests(unittest.TestCase):
    def test_two_runs_of_one_hic_tube_stay_separate(self):
        groups = cs.group_files_by_sample([
            "/reads/OG2088G-1_HICL_260420_S3_L001_R1_001.fastq.gz",
            "/reads/OG2088G-1_HICL_260420_S3_L001_R2_001.fastq.gz",
            "/reads/OG2088G-1_HICL_260605_S3_L001_R1_001.fastq.gz",
            "/reads/OG2088G-1_HICL_260605_S3_L001_R2_001.fastq.gz",
        ])
        self.assertEqual(len(groups), 2)
        self.assertEqual(
            sorted(key[2] for key in groups), ["260420", "260605"]
        )
        for files in groups.values():
            self.assertEqual(len(files["R1"]), 1)
            self.assertEqual(len(files["R2"]), 1)

    def test_lanes_of_one_run_stay_together(self):
        groups = cs.group_files_by_sample([
            "/reads/OG2088G-1_HICL_260605_S3_L001_R1_001.fastq.gz",
            "/reads/OG2088G-1_HICL_260605_S3_L002_R1_001.fastq.gz",
        ])
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(list(groups.values())[0]["R1"]), 2)

    def test_two_ilmn_dates_of_one_specimen_stay_separate(self):
        # OG768 has both 240716 and 250520 in draft-genomes.
        groups = cs.group_files_by_sample([
            "/reads/OG768.ilmn.240716.R1.fastq.gz",
            "/reads/OG768.ilmn.240716.R2.fastq.gz",
            "/reads/OG768.ilmn.250520.R1.fastq.gz",
            "/reads/OG768.ilmn.250520.R2.fastq.gz",
        ])
        self.assertEqual(len(groups), 2)
        self.assertEqual(sorted(key[2] for key in groups), ["240716", "250520"])

    def test_hic_and_ilmn_for_one_specimen_stay_separate(self):
        groups = cs.group_files_by_sample([
            "/reads/OG2088.ilmn.240716.R1.fastq.gz",
            "/reads/OG2088G-1_HICL_260605_S3_L001_R1_001.fastq.gz",
        ])
        self.assertEqual(len(groups), 2)
        self.assertEqual(sorted(key[1] for key in groups), ["hic", "ilmn"])

    def test_hifi_movies_of_one_specimen_stay_together(self):
        # A HiFi library is assembled from every movie, so movie dates must not
        # split it.
        groups = cs.group_files_by_sample([
            "/reads/OG104_m84154_241127_054344_s3.hifi_reads.bc2025.filt.fastq.gz",
            "/reads/OG104_m84154_250114_103159_s3.hifi_reads.bc2003.filt.fastq.gz",
        ])
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(list(groups.values())[0]["single"]), 2)

    def test_undated_hic_files_group_as_before(self):
        groups = cs.group_files_by_sample([
            "/reads/OG820M-1_HICL_S4_L001_R1_001.fastq.gz",
            "/reads/OG820M-1_HICL_S4_L002_R1_001.fastq.gz",
        ])
        self.assertEqual(len(groups), 1)
        self.assertEqual(list(groups)[0], ("OG820M-1_HICL", "hic", ""))


class HifiCompletionDateTests(unittest.TestCase):
    """The completion date sits at a different field index per naming style."""

    def test_single_prefix(self):
        self.assertEqual(
            cs.extract_hifi_completion_date(
                "OG785_m84154_241004_105305_s3.hifi_reads.bc2068.filt.fastq.gz"),
            "241004")

    def test_doubled_prefix(self):
        # Used to return 'm84154', which sent query_hifi_date down its except
        # ValueError branch and dated OG104's 241127 movie as 250113.
        self.assertEqual(
            cs.extract_hifi_completion_date(
                "OG104_OG104_m84154_241127_054344_s3.hifi_reads.bc2025.filt.fastq.gz"),
            "241127")

    def test_tissue_suffixed_prefix(self):
        self.assertEqual(
            cs.extract_hifi_completion_date(
                "OG62G_D_m84154_250307_141311_s4.hifi_reads.bc2010.filt.fastq.gz"),
            "250307")

    def test_sequel_movie_id(self):
        self.assertEqual(
            cs.extract_hifi_completion_date(
                "OG10G_m64497e_230209_083103.hifi_reads.bc2001--bc2001.filt.fastq.gz"),
            "230209")


class FakeCursor:
    """Answers the two queries create_samplesheet issues, from canned rows.

    `runs` maps og_id -> list of (seq_date, seq_type).
    """

    def __init__(self, runs):
        self.runs = runs
        self.result = None
        self.seen_sql = []

    def execute(self, sql, params=()):
        self.seen_sql.append(" ".join(sql.split()))
        if "FROM sequencing" in sql:
            og = params[0]
            rows = [r for r in self.runs.get(og, []) if r[1] == "PacBio"] \
                if "seq_type = 'PacBio'" in sql else list(self.runs.get(og, []))
            if len(params) == 3:
                lo, hi = params[1], params[2]
                rows = [r for r in rows if lo <= r[0] <= hi]
            rows.sort(key=lambda r: r[0], reverse=True)
            self.result = [(r[0],) for r in rows]
        elif "sample_q" in sql:
            self.result = [("Chaunax sp.", "Actinopteri", "Chaunacidae", "Lophiiformes",
                            "Chaunax")]
        else:
            self.result = []

    def fetchone(self):
        return self.result[0] if self.result else None

    def fetchall(self):
        return self.result

    def close(self):
        pass


class FakeConnection:
    def __init__(self, cursor):
        self._cursor = cursor

    def cursor(self):
        return self._cursor

    def close(self):
        pass


def run_main(tmpdir, files, runs):
    """Run create_samplesheet.main() against a fake DB and return its rows."""
    import csv
    out = Path(tmpdir) / "samplesheet.csv"
    cfg = Path(tmpdir) / "db.cfg"
    cfg.write_text("[postgres]\ndbname = x\nuser = x\npassword = x\nhost = x\nport = 5432\n")
    cursor = FakeCursor(runs)
    argv = ["create_samplesheet.py", "--output", str(out),
            "--sql-config", str(cfg), "--input-files"] + files
    with mock.patch.object(sys, "argv", argv), \
         mock.patch.object(cs.psycopg2, "connect", lambda **kw: FakeConnection(cursor),
                           create=True):
        cs.main()
    with open(out, newline="") as fh:
        return list(csv.DictReader(fh)), cursor


class HifiRunSeparationTests(unittest.TestCase):
    def test_two_runs_get_their_own_date_and_prefix(self):
        # OG104: DB runs 241127 and 250113, movies 241127 and 250114.
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            rows, _ = run_main(
                tmp,
                ["/reads/OG104_OG104_m84154_241127_054344_s3.hifi_reads.bc2025.filt.fastq.gz",
                 "/reads/OG104_OG104_m84154_250114_103159_s3.hifi_reads.bc2003.filt.fastq.gz"],
                {"OG104": [("241127", "PacBio"), ("250113", "PacBio")]})
        self.assertEqual(len(rows), 2)
        self.assertEqual({r["assembly_prefix"] for r in rows},
                         {"OG104.hifi.241127", "OG104.hifi.250113"})
        self.assertEqual({r["completion_date"] for r in rows}, {"241127", "250114"})
        self.assertEqual({r["sample"] for r in rows}, {"OG104"})

    def test_cells_of_one_run_keep_one_date(self):
        # OG2317: movies a day apart, one run. These must still be merged.
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            rows, _ = run_main(
                tmp,
                ["/reads/OG2317_m84154_260303_085443_s2.hifi_reads.bc2025.filt.fastq.gz",
                 "/reads/OG2317_m84154_260304_104941_s1.hifi_reads.bc2031.filt.fastq.gz"],
                {"OG2317": [("260303", "PacBio")]})
        self.assertEqual(len(rows), 2)
        self.assertEqual({r["assembly_prefix"] for r in rows}, {"OG2317.hifi.260303"})
        self.assertEqual({r["date"] for r in rows}, {"260303"})

    def test_unreadable_completion_date_is_unknown_not_the_latest_run(self):
        # An unparseable date used to fall back to the sample's most recent
        # PacBio run, which is how OG104's 241127 movie was labelled 250113.
        cursor = FakeCursor({"OG104": [("241127", "PacBio"), ("250113", "PacBio")]})
        self.assertEqual(cs.query_hifi_date(cursor, "OG104", "m84154"), "unknown")
        self.assertEqual(cursor.seen_sql, [])

    def test_hifi_is_never_dated_from_another_technology(self):
        # OG16's PacBio run is outside the window; a HiC run inside it must not
        # be borrowed.
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            rows, cursor = run_main(
                tmp,
                ["/reads/OG16_m84154_240417_083812_s1.hifi_reads.bc2079.filt.fastq.gz"],
                {"OG16": [("240405", "HiC"), ("240427", "PacBio")]})
        self.assertEqual(rows[0]["date"], "unknown")
        self.assertTrue(any("seq_type = 'PacBio'" in sql for sql in cursor.seen_sql))


if __name__ == "__main__":
    unittest.main()
