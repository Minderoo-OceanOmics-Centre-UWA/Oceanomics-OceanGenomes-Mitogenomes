"""Unit tests for the [mgcode=] tag the gene/CDS extractors write.

Both extractors used to hardcode `[mgcode=2]` into every FASTA header they
emitted, while the genome FASTA, the .tbl and table2asn all carried the
sample's real translation table. A code-4 coral therefore shipped an ENA
candidate package whose gene and protein FASTAs contradicted its own genome
record. These tests pin the header to the code the caller passes, and pin the
refusal to run without one.

Both scripts import Biopython at module scope, so skip cleanly when it is
absent (same pattern as test_mitos_to_emma.py).
"""

import subprocess
import sys
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"

try:
    import Bio  # noqa: F401
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


def run(script, args, cwd):
    """Run a bin/ script the way Nextflow does: siblings flat on PATH."""
    env = {"PATH": f"{BIN}:/usr/bin:/bin", "PYTHONPATH": str(BIN)}
    return subprocess.run(
        [sys.executable, str(BIN / script), *args],
        cwd=cwd, env=env, capture_output=True, text=True,
    )


GENOME = "ATG" + "AAA" * 20 + "TAA"


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class ExtractGenesGffTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        (self.root / "asm.fa").write_text(f">OG1.ilmn.240101.getorg1770\n{GENOME}\n")
        (self.root / "asm.gff").write_text(textwrap.dedent("""\
            ##gff-version 3
            ##organism Acropora millepora
            asm\ttest\tgene\t1\t66\t.\t+\t.\tName=CO1;Product=cytochrome c oxidase subunit I
            asm\ttest\tCDS\t1\t66\t.\t+\t.\tName=CO1;Product=cytochrome c oxidase subunit I
            """))

    def extract(self, code):
        outdir = self.root / f"out{code}"
        res = run("extract_genes_gff.py", [
            "--fasta", "asm.fa", "--gff", "asm.gff",
            "--outdir", str(outdir), "--assembly", "OG1.ilmn.240101.getorg1770",
            "--genetic-code", str(code),
        ], cwd=self.root)
        return res, outdir

    def test_coral_header_carries_code_4_not_2(self):
        res, outdir = self.extract(4)
        self.assertEqual(res.returncode, 0, res.stderr)
        text = (outdir / "genes" / "OG1.ilmn.240101.getorg1770.genes.fa").read_text()
        self.assertIn("[mgcode=4]", text)
        self.assertNotIn("[mgcode=2]", text)

    def test_vertebrate_header_still_carries_code_2(self):
        res, outdir = self.extract(2)
        self.assertEqual(res.returncode, 0, res.stderr)
        text = (outdir / "genes" / "OG1.ilmn.240101.getorg1770.genes.fa").read_text()
        self.assertIn("[mgcode=2]", text)

    def test_per_cds_single_matches_the_concatenated_header(self):
        _res, outdir = self.extract(9)
        single = (outdir / "genes" / "CO1.OG1.ilmn.240101.getorg1770.fa").read_text()
        self.assertIn("[mgcode=9]", single)

    def test_missing_genetic_code_is_refused(self):
        res = run("extract_genes_gff.py", [
            "--fasta", "asm.fa", "--gff", "asm.gff",
            "--outdir", str(self.root / "none"), "--assembly", "OG1",
        ], cwd=self.root)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("--genetic-code", res.stderr)

    def test_unsupported_genetic_code_is_refused(self):
        res, _outdir = self.extract(7)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("not a supported", res.stdout + res.stderr)


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class ExtractCdsFromTblTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        (self.root / "asm.fa").write_text(
            f">OG1.ilmn.240101.getorg1770 [organism=Acropora millepora]\n{GENOME}\n")
        (self.root / "asm.tbl").write_text(
            "gene\tCO1\n1\t66\tCDS\n\t\t\tproduct\tcytochrome c oxidase subunit I\n")

    def extract(self, code):
        outdir = self.root / f"out{code}"
        res = run("extract_cds_from_tbl.py", [
            "--fasta", "asm.fa", "--tbl", "asm.tbl",
            "--outdir", str(outdir), "--assembly", "OG1.ilmn.240101.getorg1770",
            "--genetic-code", str(code),
        ], cwd=self.root)
        return res, outdir

    def test_coral_cds_header_carries_code_4(self):
        res, outdir = self.extract(4)
        self.assertEqual(res.returncode, 0, res.stderr)
        text = (outdir / "cds" / "CO1.OG1.ilmn.240101.getorg1770.fa").read_text()
        self.assertIn("[mgcode=4]", text)
        self.assertNotIn("[mgcode=2]", text)

    def test_missing_genetic_code_is_refused(self):
        res = run("extract_cds_from_tbl.py", [
            "--fasta", "asm.fa", "--tbl", "asm.tbl",
            "--outdir", str(self.root / "none"), "--assembly", "OG1",
        ], cwd=self.root)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("--genetic-code", res.stderr)


if __name__ == "__main__":
    unittest.main()
