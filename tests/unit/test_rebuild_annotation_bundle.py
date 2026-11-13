"""Unit tests for bin/rebuild_annotation_bundle.py.

The --blast_precomputed_annotations path rebuilds each precomputed bundle's cds/ and
proteins/ from its GFF + FASTA (a glob cannot carry those directories, and without them
ANNOTATION_STATS skips the PCG ORF check) and writes the CO1 / 12S / 16S BLAST inputs.
These tests pin:

  * real MITOS and EMMA bundles reproduce the cds/, proteins/ and BLAST inputs their
    annotators wrote, so a re-run checks and BLASTs exactly what the original run did;
  * a protein that has to be re-translated is byte-identical to what MITOS publishes;
  * the filename/header contract that getAnnotationName(), the BLAST/LCA names and
    annotation_stats.py rely on;
  * origin-spanning and multi-exon features are reconstructed, not truncated;
  * a missing region is recorded, not invented;
  * a curated GFF beats a stale cds/, and the mismatch is flagged.

Biopython is imported at module scope by the script, so skip cleanly when it is absent
(same pattern as test_extract_genes_mgcode.py).
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"
FIXTURES = Path(__file__).resolve().parent / "fixtures" / "annotation_bundles"

MITOS_STEM = "OG2905.ilmn.260925.getorg1770.mitos2110"   # Asteroidea, genetic code 9
EMMA_STEM = "OG3148.ilmn.260925.getorg1770.emma102"      # Glaucosoma hebraicum, code 2

try:
    import Bio  # noqa: F401
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


def run(files, workdir, code=None):
    """Run the script the way Nextflow does: siblings flat on PATH."""
    env = {"PATH": f"{BIN}:/usr/bin:/bin", "PYTHONPATH": str(BIN)}
    cmd = [sys.executable, str(BIN / "rebuild_annotation_bundle.py"),
           "--outdir", str(workdir / "annotation"), "--regions-dir", str(workdir / "regions")]
    if code is not None:
        cmd += ["--genetic-code", str(code)]
    return subprocess.run(cmd + [str(f) for f in files], cwd=workdir, env=env,
                          capture_output=True, text=True)


def read_fasta(path):
    lines = Path(path).read_text().splitlines()
    return lines[0], "".join(lines[1:])


def status(workdir, stem):
    with open(workdir / "regions" / f"{stem}.bundle_rebuild.tsv") as fh:
        return {row["gene"]: row for row in csv.DictReader(fh, delimiter="\t")}


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class RealBundleTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = Path(self._tmp.name)

    def stage(self, bundle, symlink=True):
        """Stage the bundle's top-level files as Nextflow does (symlinks, so cds/ and
        proteins/ are found beside the real GFF), or as copies to hide them."""
        stage = self.tmp / "bundle"
        stage.mkdir()
        for f in (FIXTURES / bundle).iterdir():
            if f.is_file():
                (os.symlink if symlink else shutil.copy)(f, stage / f.name)
        return sorted(stage.iterdir())

    def assert_same_seq(self, ours, theirs, label):
        self.assertEqual(read_fasta(ours)[1], read_fasta(theirs)[1], label)

    def check_bundle(self, bundle, stem, code, cds_prefix, expected_check):
        res = run(self.stage(bundle), self.tmp, code)
        self.assertEqual(res.returncode, 0, res.stderr)
        rows = status(self.tmp, stem)
        out = self.tmp / "annotation"
        # Top-level files carried over, so FORMAT_FILES / ANNOTATION_STATS still find them.
        self.assertTrue((out / f"{stem}.fa").is_file())
        self.assertTrue((out / f"{stem}.gff").is_file())
        pcgs = [g for g, r in rows.items() if r["kind"] == "CDS"]
        self.assertEqual(len(pcgs), 13)
        self.assertEqual(len(rows), 15)
        for gene, row in rows.items():
            self.assertEqual(row["status"], "found", gene)
            self.assertEqual(row["cds_check"], expected_check, gene)
            self.assert_same_seq(out / "cds" / f"{gene}.{stem}.fa",
                                 FIXTURES / bundle / "cds" / f"{cds_prefix}{gene}.{stem}.fa", gene)
        for gene in pcgs:
            self.assert_same_seq(out / "proteins" / f"MT-{gene}.{stem}.fa",
                                 FIXTURES / bundle / "proteins" / f"MT-{gene}.{stem}.fa", gene)
        for region in ("CO1", "RNR1", "RNR2"):
            header, _seq = read_fasta(self.tmp / "regions" / f"{region}.{stem}.fa")
            self.assertEqual(header, f">{stem}")
            self.assert_same_seq(self.tmp / "regions" / f"{region}.{stem}.fa",
                                 FIXTURES / bundle / "cds" / f"{cds_prefix}{region}.{stem}.fa", region)

    def test_mitos_bundle_matches_published_files(self):
        self.check_bundle("mitos", MITOS_STEM, 9, "", "match")

    def test_emma_bundle_reuses_its_trimmed_cds_and_proteins(self):
        # EMMA's cds/ drops the stop codon and the last base or two of each rRNA.
        self.check_bundle("emma", EMMA_STEM, 2, "MT-", "trimmed")

    def test_retranslated_proteins_match_mitos(self):
        # Copies, not symlinks: no cds/ or proteins/ to reuse, so every protein is
        # translated from the GFF span and must equal what MITOS published.
        res = run(self.stage("mitos", symlink=False), self.tmp, 9)
        self.assertEqual(res.returncode, 0, res.stderr)
        rows = status(self.tmp, MITOS_STEM)
        for gene, row in rows.items():
            self.assertEqual(row["cds_check"], "no_cds", gene)
            if row["kind"] == "CDS":
                self.assertIn("protein_translated", row["notes"])
                self.assert_same_seq(
                    self.tmp / "annotation" / "proteins" / f"MT-{gene}.{MITOS_STEM}.fa",
                    FIXTURES / "mitos" / "proteins" / f"MT-{gene}.{MITOS_STEM}.fa", gene)

    def test_retranslation_needs_a_genetic_code(self):
        res = run(self.stage("mitos", symlink=False), self.tmp)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("--genetic-code", res.stderr)

    def test_two_gffs_in_one_bundle_is_an_error(self):
        files = self.stage("mitos")
        extra = self.tmp / "bundle" / f"{MITOS_STEM}_mancurated.gff"
        shutil.copy(FIXTURES / "mitos" / f"{MITOS_STEM}.gff", extra)
        res = run(files + [extra], self.tmp, 9)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("more than one GFF", res.stderr)


# 60 bp genome: positions 1-based. Distinct blocks so slices are easy to read.
GENOME = "AAAAACCCCC" "GGGGGTTTTT" "ACGTACGTAC" "TTTTTGGGGG" "CCCCCAAAAA" "GATCGATCGA"


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython not installed")
class SyntheticTests(unittest.TestCase):
    STEM = "OG1.ilmn.260101.getorg1770.mitos2110"

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp = Path(self._tmp.name)
        (self.tmp / f"{self.STEM}.fa").write_text(f">OG1.ilmn.260101.getorg1770\n{GENOME}\n")

    def rebuild(self, gff_rows, cds=None):
        gff = self.tmp / f"{self.STEM}.gff"
        gff.write_text("##gff-version 3\n" + textwrap.dedent(gff_rows))
        if cds:
            (self.tmp / "cds").mkdir()
            for gene, seq in cds.items():
                (self.tmp / "cds" / f"{gene}.{self.STEM}.fa").write_text(f">{self.STEM}\n{seq}\n")
        res = run([self.tmp / f"{self.STEM}.fa", gff], self.tmp, code=2)
        self.assertEqual(res.returncode, 0, res.stderr)
        return status(self.tmp, self.STEM), res

    def region(self, gene):
        return read_fasta(self.tmp / "regions" / f"{gene}.{self.STEM}.fa")[1]

    def test_origin_spanning_feature_joins_tail_and_head(self):
        self.rebuild("""\
            c\tmitos\tgene\t56\t60\t.\t+\t.\tID=gene-CO1;Name=MT-CO1;Note=origin-spanning;wrapped_end=5
            c\tmitos\tCDS\t56\t60\t.\t+\t.\tID=CO1;Parent=gene-CO1;Name=MT-CO1;Note=origin-spanning;wrapped_end=5
            """)
        self.assertEqual(self.region("CO1"), GENOME[55:] + GENOME[:5])

    def test_multi_exon_gene_is_spliced_in_transcript_order(self):
        rows, _res = self.rebuild("""\
            c\tmitos\tgene\t1\t30\t.\t+\t.\tID=gene-CO1;Name=MT-CO1;Note=multi-exon
            c\tmitos\tCDS\t1\t10\t.\t+\t.\tID=CO1.1;Parent=gene-CO1;Name=MT-CO1
            c\tmitos\tCDS\t21\t30\t.\t+\t.\tID=CO1.2;Parent=gene-CO1;Name=MT-CO1
            """)
        self.assertEqual(self.region("CO1"), GENOME[0:10] + GENOME[20:30])
        self.assertIn("exons=2", rows["CO1"]["notes"])

    def test_minus_strand_rrna_is_reverse_complemented(self):
        self.rebuild("""\
            c\tmitos\trRNA\t11\t20\t.\t-\t.\tID=RNR2;Parent=gene-RNR2;Name=MT-RNR2
            """)
        self.assertEqual(self.region("RNR2"), "AAAAACCCCC")  # revcomp of GGGGGTTTTT

    def test_missing_region_is_recorded_not_written(self):
        rows, res = self.rebuild("""\
            c\tmitos\tCDS\t1\t10\t.\t+\t.\tID=CO1;Parent=gene-CO1;Name=MT-CO1
            c\tmitos\trRNA\t21\t30\t.\t+\t.\tID=RNR2;Parent=gene-RNR2;Name=MT-RNR2
            """)
        self.assertFalse((self.tmp / "regions" / f"RNR1.{self.STEM}.fa").exists())
        self.assertEqual(rows["RNR1"]["status"], "missing")
        self.assertIn("2/3 BLAST regions", res.stdout)

    def test_mitos_alias_names_are_recognised_and_trnas_skipped(self):
        rows, _res = self.rebuild("""\
            c\tx\tCDS\t1\t9\t.\t+\t.\tName=cox1
            c\tx\trRNA\t11\t20\t.\t+\t.\tName=rrnS
            c\tx\trRNA\t21\t30\t.\t+\t.\tName=16S
            c\tx\ttRNA\t31\t40\t.\t+\t.\tName=trnF(gaa)
            """)
        self.assertEqual(self.region("CO1"), GENOME[0:9])
        self.assertEqual(self.region("RNR1"), GENOME[10:20])
        self.assertEqual(self.region("RNR2"), GENOME[20:30])
        self.assertNotIn("TF", rows)
        # Only protein-coding genes get a protein.
        self.assertTrue((self.tmp / "annotation" / "proteins" / f"MT-CO1.{self.STEM}.fa").is_file())
        self.assertFalse((self.tmp / "annotation" / "proteins" / f"MT-RNR1.{self.STEM}.fa").exists())

    def test_curated_gff_beats_stale_cds(self):
        rows, res = self.rebuild("""\
            c\tmitos\tCDS\t1\t21\t.\t+\t.\tID=CO1;Parent=gene-CO1;Name=MT-CO1
            """, cds={"CO1": GENOME[0:30]})   # cds/ from before the curation shortened CO1
        self.assertEqual(self.region("CO1"), GENOME[0:21])
        self.assertEqual(read_fasta(self.tmp / "annotation" / "cds" / f"CO1.{self.STEM}.fa")[1],
                         GENOME[0:21])
        self.assertEqual(rows["CO1"]["cds_check"], "stale")
        self.assertIn("protein_translated", rows["CO1"]["notes"])
        self.assertIn("WARNING", res.stderr)

    def test_no_cds_dir_uses_the_gff(self):
        rows, _res = self.rebuild("""\
            c\tmitos\tCDS\t1\t21\t.\t+\t.\tID=CO1;Parent=gene-CO1;Name=MT-CO1
            """)
        self.assertEqual(self.region("CO1"), GENOME[0:21])
        self.assertEqual(rows["CO1"]["cds_check"], "no_cds")


if __name__ == "__main__":
    unittest.main()
