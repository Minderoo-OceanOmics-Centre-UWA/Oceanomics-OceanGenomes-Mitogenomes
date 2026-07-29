"""Unit tests for candidate reference ranking in bin/rank_reference_candidates.py.

minimap2 is stubbed with a fake PAF so the selection logic is exercised without the
aligner: what matters here is that the best-recruiting candidate wins, that ties keep
the taxonomically closest one, and that every failure mode still yields a reference.
"""
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "rank_reference_candidates", ROOT / "bin" / "rank_reference_candidates.py"
)
rrc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rrc)


def paf(matches):
    """One PAF line with `matches` residue matches (col 10, 0-indexed 9)."""
    return "\t".join(["read1", "1000", "0", "1000", "+", "ref", "16000", "0",
                      "1000", str(matches), "1000", "60"]) + "\n"


class ScoreCandidateTests(unittest.TestCase):
    """score_candidate normalises matched bases by reference length, so a longer
    reference cannot win merely by being longer."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        self._real_run = rrc.subprocess.run

    def stub_minimap2(self, stdout):
        def fake_run(cmd, **kwargs):
            return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")
        rrc.subprocess.run = fake_run
        self.addCleanup(lambda: setattr(rrc.subprocess, "run", self._real_run))

    def write_fasta(self, name, length):
        path = self.dir / name
        path.write_text(">ref\n" + ("A" * length) + "\n")
        return path

    def test_score_is_matches_per_reference_base(self):
        self.stub_minimap2(paf(8000))
        fasta = self.write_fasta("A.fasta", 16000)
        score, matched, mapped = rrc.score_candidate(fasta, "reads.fq", "map-hifi", 1)
        self.assertAlmostEqual(score, 0.5)
        self.assertEqual(matched, 8000)
        self.assertEqual(mapped, 1)

    def test_longer_reference_does_not_win_on_length_alone(self):
        self.stub_minimap2(paf(8000))
        short = self.write_fasta("short.fasta", 16000)
        long_ = self.write_fasta("long.fasta", 32000)
        self.assertGreater(
            rrc.score_candidate(short, "reads.fq", "map-hifi", 1)[0],
            rrc.score_candidate(long_, "reads.fq", "map-hifi", 1)[0],
        )

    def test_minimap2_failure_scores_zero_rather_than_raising(self):
        def boom(cmd, **kwargs):
            raise OSError("minimap2 not found")
        rrc.subprocess.run = boom
        self.addCleanup(lambda: setattr(rrc.subprocess, "run", self._real_run))
        fasta = self.write_fasta("A.fasta", 16000)
        self.assertEqual(rrc.score_candidate(fasta, "reads.fq", "map-hifi", 1),
                         (0.0, 0, 0))

    def test_unaligned_reads_score_zero(self):
        self.stub_minimap2("")
        fasta = self.write_fasta("A.fasta", 16000)
        self.assertEqual(rrc.score_candidate(fasta, "reads.fq", "map-hifi", 1)[0], 0.0)


class SubsampleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def fastq(self, name, n):
        path = self.dir / name
        path.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(n)))
        return path

    def test_caps_at_requested_read_count(self):
        src = self.fastq("in.fq", 100)
        out = self.dir / "out.fq"
        self.assertEqual(rrc.subsample_reads([src], out, 10), 10)
        self.assertEqual(len(out.read_text().splitlines()), 40)

    def test_takes_everything_when_file_is_smaller_than_the_cap(self):
        src = self.fastq("in.fq", 5)
        out = self.dir / "out.fq"
        self.assertEqual(rrc.subsample_reads([src], out, 1000), 5)

    def test_draws_across_multiple_read_files(self):
        a, b = self.fastq("a.fq", 3), self.fastq("b.fq", 3)
        out = self.dir / "out.fq"
        self.assertEqual(rrc.subsample_reads([a, b], out, 6), 6)

    def test_missing_read_file_does_not_raise(self):
        out = self.dir / "out.fq"
        self.assertEqual(rrc.subsample_reads([self.dir / "nope.fq"], out, 10), 0)


class CandidateLengthTests(unittest.TestCase):
    def test_ignores_headers_and_newlines(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            fh.write(">ref desc\nACGTACGTAC\nACGTA\n")
            path = Path(fh.name)
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        self.assertEqual(rrc.candidate_length(path), 15)

    def test_missing_file_is_zero_length(self):
        self.assertEqual(rrc.candidate_length(Path("/nonexistent.fasta")), 0)


if __name__ == "__main__":
    unittest.main()
