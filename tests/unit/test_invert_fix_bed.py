"""Reference-free ORF repair, available to every invertebrate.

The snap itself is not new; it was locked inside coral_fix_bed.py, so in
NOVA_260909_LA all 17 assemblies carrying a genuinely broken CDS -- 13 Porifera,
the rest Malacostraca and Asteroidea, not one of them a coral -- had no repair
path at all. Five Demospongiae shared one systematic MITOS boundary miscall.

What matters in these tests is the refusals as much as the repairs: this script
reads no reference, so it can never import a gene the assembly does not carry,
and where it declines it must say why rather than leaving a silent no-op.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))

import orf_snap
from orf_utils import classify_cds

SCRIPT = BIN / "invert_fix_bed.py"


def run_fixer(bed_rows, genome, code, lineage=None):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "in.bed").write_text(
            "\n".join("\t".join(str(c) for c in r) for r in bed_rows) + "\n")
        (tmp / "g.fa").write_text(">c\n" + genome + "\n")
        cmd = [sys.executable, str(SCRIPT), "--bed", str(tmp / "in.bed"),
               "--genome", str(tmp / "g.fa"), "--code", str(code),
               "--out-bed", str(tmp / "out.bed"), "--status", str(tmp / "s.tsv")]
        if lineage:
            cmd += ["--lineage", lineage]
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"})
        return proc, (tmp / "out.bed").read_text(), (tmp / "s.tsv").read_text()


class SharedMachineryTests(unittest.TestCase):
    def test_nad5_is_snappable(self):
        """nad5 was excluded as a module constant because the CORAL reference pass
        owns it. For every non-coral invertebrate nothing else touches nad5, so
        the exclusion left it permanently unrepairable outside Anthozoa."""
        self.assertIn("nad5", orf_snap.PCG_BARE)

    def test_coral_fixer_uses_the_same_implementation(self):
        import coral_fix_bed
        self.assertIs(coral_fix_bed.snap_pcg_rows, orf_snap.snap_pcg_rows)

    def test_coral_fixer_still_excludes_nad5_by_passing_skip(self):
        """The exclusion is the caller's to declare, not a property of nad5."""
        text = (BIN / "coral_fix_bed.py").read_text()
        self.assertIn('skip=({"nad5", "cox1"}', text)


class BehaviourTests(unittest.TestCase):
    def test_always_exits_zero(self):
        """A repair attempt, not a gate. The QC step decides publishability."""
        proc, _bed, _status = run_fixer(
            [["c", 0, 60, "cox1", 0, "+"]], "ATG" + "GGG" * 19 + "TAA", 5)
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_clean_annotation_is_left_untouched(self):
        genome = "ATG" + "GGG" * 19 + "TAA" + "TTTT"
        rows = [["c", 0, 63, "cox1", 0, "+"]]
        _proc, out_bed, status = run_fixer(rows, genome, 5)
        self.assertIn("SKIPPED", status.splitlines()[0])
        self.assertIn("cox1", out_bed)

    def test_status_names_every_gene_it_could_not_repair(self):
        """A boundary left alone is a fact the QC step needs, not an absence."""
        genome = "TTT" + "GGG" * 19 + "AGA" + "CCCCCCCCC"
        rows = [["c", 0, 63, "nad4l", 0, "+"]]
        _proc, _bed, status = run_fixer(rows, genome, 4)
        self.assertRegex(status, r"(FAILED|PARTIAL)")
        self.assertIn("nad4l", status)

    def test_output_bed_is_always_written(self):
        _proc, out_bed, _status = run_fixer(
            [["c", 0, 60, "nad1", 0, "+"]], "A" * 200, 5)
        self.assertTrue(out_bed.strip(), "a BED must be emitted even when nothing changed")

    def test_documented_alternative_start_is_not_treated_as_broken(self):
        """A decapod CO1 opening on CGA is correct as annotated; 'repairing' it
        would move a boundary that was already right."""
        genome = "CGA" + "GGG" * 19 + "TAA" + "TTTT"
        rows = [["c", 0, 63, "cox1", 0, "+"]]
        _proc, _bed, status = run_fixer(rows, genome, 5, lineage="Malacostraca")
        self.assertIn("SKIPPED", status.splitlines()[0])
        self.assertNotIn("PCG-fix cox1", status)

    def test_same_codon_without_the_lineage_is_still_examined(self):
        genome = "CGA" + "GGG" * 19 + "TAA" + "TTTT"
        rows = [["c", 0, 63, "cox1", 0, "+"]]
        _proc, _bed, status = run_fixer(rows, genome, 5)
        self.assertNotIn("SKIPPED", status.splitlines()[0])


class LinearTopologyTests(unittest.TestCase):
    """A linear contig's ends are assembly breaks. A gene cut by the break must not
    be 'repaired' with a stop codon borrowed from the other end of the contig."""

    # The only in-frame stop for the gene at the 3' end sits at the very start of
    # the contig, reachable only by wrapping the origin.
    GENE = "ATG" + "GCT" * 49            # 150 nt, no stop
    GENOME = "TAA" + "CCC" * 40 + GENE
    ROW = ["c", len(GENOME) - 150, len(GENOME), "atp8", 0, "+"]

    def test_circular_snap_may_wrap(self):
        fixes, _actions, _unrep = orf_snap.snap_pcg_rows(
            [list(self.ROW)], self.GENOME, 5, circular=True)
        self.assertEqual(len(fixes), 1)
        s0, e1 = fixes[0]
        self.assertLess(e1, s0, "circular snap should take the stop across the origin")

    def test_linear_snap_never_wraps(self):
        fixes, _actions, unrep = orf_snap.snap_pcg_rows(
            [list(self.ROW)], self.GENOME, 5, circular=False)
        for s0, e1 in fixes.values():
            self.assertLess(s0, e1, "a linear contig has no origin to cross")
        self.assertIn("atp8", unrep)

    def test_linear_flag_reaches_the_snap(self):
        proc, out_bed, status = run_fixer([self.ROW], self.GENOME, 5)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "in.bed").write_text("\t".join(str(c) for c in self.ROW) + "\n")
            (tmp / "g.fa").write_text(">c\n" + self.GENOME + "\n")
            proc = subprocess.run(
                [sys.executable, str(SCRIPT), "--bed", str(tmp / "in.bed"),
                 "--genome", str(tmp / "g.fa"), "--code", "5", "--linear",
                 "--out-bed", str(tmp / "out.bed"), "--status", str(tmp / "s.tsv")],
                capture_output=True, text=True,
                env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            f = (tmp / "out.bed").read_text().split("\t")
            self.assertLess(int(f[1]), int(f[2]))


class StatusFormatTests(unittest.TestCase):
    def test_status_has_a_verdict_and_a_per_gene_table(self):
        genome = "TTT" + "GGG" * 19 + "AGA" + "CCCCCCCCC"
        _proc, _bed, status = run_fixer([["c", 0, 63, "nad4l", 0, "+"]], genome, 4)
        lines = status.splitlines()
        self.assertIn(lines[0].split("\t")[0],
                      {"FIXED", "PARTIAL", "SKIPPED", "FAILED"})
        self.assertTrue(any(l.startswith("#gene") for l in lines))


if __name__ == "__main__":
    unittest.main()
