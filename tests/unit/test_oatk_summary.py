"""Unit tests for the Oatk (reference-free HiFi fallback) branch of the assembly
summary: run classification, prefix inference, GFA circularity, and end-to-end
parsing of an Oatk run into a complete row."""
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "mitogenome_assembly_summary", ROOT / "bin" / "mitogenome_assembly_summary.py"
)
mas = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mas
SPEC.loader.exec_module(mas)

PFX = "OG62.hifi.260227.v10oatk"


class OatkClassificationTests(unittest.TestCase):
    def test_classify_oatk_files(self):
        self.assertEqual(mas.classify_file(Path(f"x/{PFX}.fasta")), "Oatk")
        self.assertEqual(mas.classify_file(Path(f"x/{PFX}.mito.ctg.fasta")), "Oatk")
        self.assertEqual(mas.classify_file(Path(f"x/{PFX}.mito.gfa")), "Oatk")
        self.assertEqual(mas.classify_file(Path(f"x/{PFX}.oatk.log")), "Oatk")

    def test_oatk_not_confused_with_mitohifi(self):
        # A real mitohifi file must still classify as MitoHiFi.
        self.assertEqual(mas.classify_file(Path("x/OG1.hifi.260101.v323mitohifi.fasta")), "MitoHiFi")

    def test_infer_prefix_consistent_across_oatk_files(self):
        for name in (f"{PFX}.fasta", f"{PFX}.mito.ctg.fasta", f"{PFX}.mito.gfa", f"{PFX}.oatk.log"):
            self.assertEqual(mas.infer_prefix(Path(name), "Oatk"), PFX, name)


class OatkGfaCircularityTests(unittest.TestCase):
    def _gfa(self, text):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        p = Path(tmp.name) / f"{PFX}.mito.gfa"
        p.write_text(text)
        return [p]

    def test_self_link_is_circular(self):
        files = self._gfa("S\tu0\tACGT\tLN:i:16465\nL\tu0\t+\tu0\t+\t985M\n")
        self.assertEqual(mas.oatk_circular_from_gfa(files), "true")

    def test_no_self_link_is_linear(self):
        files = self._gfa("S\tu0\tACGT\nS\tu1\tACGT\nL\tu0\t+\tu1\t+\t0M\n")
        self.assertEqual(mas.oatk_circular_from_gfa(files), "false")

    def test_no_gfa_unknown(self):
        self.assertEqual(mas.oatk_circular_from_gfa([Path("x/other.txt")]), "")


class OatkParseRunTests(unittest.TestCase):
    def build_run(self, circular=True, genes="37", cds="13"):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        run = root / "OG62" / PFX / "mtdna"
        ann = root / "OG62" / PFX / "annotation"
        run.mkdir(parents=True); ann.mkdir(parents=True)
        (run / f"{PFX}.fasta").write_text(">c\n" + "ACGT" * 4116 + "A\n")  # 16465 bp
        link = "L\tu0\t+\tu0\t+\t985M\n" if circular else "L\tu0\t+\tu1\t+\t0M\n"
        (run / f"{PFX}.mito.gfa").write_text("S\tu0\tACGT\n" + link)
        (run / f"{PFX}.oatk.log").write_text("oatk log\n")
        (ann / f"{PFX}.annotation_stats.csv").write_text(
            f"sample,num_genes,num_cds,missing_genes,frameshift_flag\n{PFX},{genes},{cds},no,false\n"
        )
        return root

    def summarise(self, root):
        thr = mas.Thresholds(20, 1.0, 10000, 25000, 37, 13)
        rows = mas.build_summary([root], thr)
        return [r for r in rows if r["assembler"] == "Oatk"][0]

    def test_complete_circular_oatk_run(self):
        row = self.summarise(self.build_run(circular=True))
        self.assertEqual(row["status"], "complete")
        self.assertEqual(row["circularised"], "true")
        self.assertEqual(row["final_length_bp"], "16465")
        self.assertEqual(row["num_final_contigs"], "1")
        self.assertEqual(row["num_genes"], "37")
        self.assertEqual(row["reference_species"], "")  # reference-free

    def test_linear_oatk_run_flagged(self):
        row = self.summarise(self.build_run(circular=False))
        self.assertEqual(row["status"], "manual_review")
        self.assertIn("not_circularised", row["manual_review_reason"])


if __name__ == "__main__":
    unittest.main()
