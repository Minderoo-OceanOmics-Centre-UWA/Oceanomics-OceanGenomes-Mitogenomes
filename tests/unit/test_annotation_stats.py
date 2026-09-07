"""Unit tests for bin/annotation_stats.py -- the completeness gate.

Covers process_gff(): missing_genes / order_correct / passed and the bounded
tRNA-tolerance branch (trna_advisory). Pure stdlib -- the Biopython import lives
inside process_protein_lengths(), which these tests do not call -- so they always
run.
"""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, mito_gene_order, ...) the way
# Nextflow stages them: flat on PATH. Mirror that for the file-path loads below.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "annotation_stats", ROOT / "bin" / "annotation_stats.py")
stats = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stats)

REF = stats.REF_GENES


def _feature_type(gene):
    if gene.startswith("RNR"):
        return "rRNA"
    if gene.startswith("T"):
        return "tRNA"
    return "CDS"


def gff_for(genes, region_len=16500):
    """Minimal EMMA-style GFF: one `gene` line + one matching feature line per
    name, at ascending coordinates in the given order."""
    lines = ["##gff-version 3", f"##sequence-region chr 1 {region_len}"]
    pos = 100
    for g in genes:
        end = pos + 50
        gid = f"gene-{g}"
        lines.append(f"chr\tEmma\tgene\t{pos}\t{end}\t.\t+\t.\tID={gid};Name=MT-{g}")
        lines.append(
            f"chr\tEmma\t{_feature_type(g)}\t{pos}\t{end}\t.\t+\t.\t"
            f"ID=feat-{g};Parent={gid};Name=MT-{g}")
        pos += 100
    return "\n".join(lines) + "\n"


class ProcessGffTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _run(self, genes, trna_tolerance=2, class_name=""):
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text(gff_for(genes))
        return stats.process_gff(str(p), p.stem, class_name, trna_tolerance)

    def test_complete_in_order_passes(self):
        s = self._run(REF)
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["missing_genes"], "no")
        self.assertEqual(s["trna_advisory"], "no")
        self.assertEqual(s["order_correct"], "yes")

    def test_one_trna_missing_tolerated(self):
        genes = [g for g in REF if g != "TP"]
        s = self._run(genes, trna_tolerance=2)
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["missing_genes"], "TP")      # stays truthful
        self.assertEqual(s["trna_advisory"], "TP")

    def test_three_trna_missing_not_tolerated_at_2(self):
        genes = [g for g in REF if g not in ("TW", "TA", "TN")]
        s = self._run(genes, trna_tolerance=2)
        self.assertEqual(s["passed"], "no")
        self.assertEqual(s["trna_advisory"], "no")
        self.assertEqual(s["missing_genes"], "TW;TA;TN")

    def test_three_trna_missing_tolerated_at_3(self):
        genes = [g for g in REF if g not in ("TW", "TA", "TN")]
        s = self._run(genes, trna_tolerance=3)
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["trna_advisory"], "TW;TA;TN")

    def test_missing_pcg_never_tolerated(self):
        genes = [g for g in REF if g != "ND4L"]
        s = self._run(genes, trna_tolerance=3)
        self.assertEqual(s["passed"], "no")
        self.assertEqual(s["trna_advisory"], "no")

    def test_missing_rrna_never_tolerated(self):
        genes = [g for g in REF if g != "RNR2"]
        s = self._run(genes, trna_tolerance=3)
        self.assertEqual(s["passed"], "no")

    def test_trna_missing_but_order_broken_not_tolerated(self):
        genes = [g for g in REF if g != "TP"]
        genes[5], genes[6] = genes[6], genes[5]   # scramble two genes
        s = self._run(genes, trna_tolerance=2)
        self.assertEqual(s["order_correct"], "no")
        self.assertEqual(s["passed"], "no")
        self.assertEqual(s["trna_advisory"], "no")

    def test_tolerance_zero_requires_complete(self):
        genes = [g for g in REF if g != "TP"]
        s = self._run(genes, trna_tolerance=0)
        self.assertEqual(s["passed"], "no")

    def test_cnidarian_branch_unaffected(self):
        # Core (13 PCG + 2 rRNA) present, tRNAs absent -> cnidarian passes,
        # trna_advisory not applicable.
        core = [g for g in REF if not g.startswith("T")]
        s = self._run(core, class_name="Anthozoa")
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["trna_advisory"], "no")
        self.assertEqual(s["order_correct"], "NA")


class CompletenessProfileTests(unittest.TestCase):
    """Which profile an assembly is judged under.

    Previously keyed on a hardcoded cnidarian class list, so a code-9 echinoderm
    or code-5 mollusc was judged against the vertebrate 37-gene set and gene
    order, and failed for being what it is. The resolved genetic code decides it
    now; the class stays as the fallback for callers that have no code.
    """

    def test_code_2_is_the_vertebrate_profile(self):
        self.assertEqual(stats.completeness_profile(2), "vertebrate")

    def test_every_other_code_is_the_core_profile(self):
        for code in (4, 5, 9, 13, 14, 21, 24, 33):
            with self.subTest(code=code):
                self.assertEqual(stats.completeness_profile(code), "core")

    def test_code_beats_a_disagreeing_class(self):
        # A mislabelled class must not drag a code-4 sample onto the vertebrate
        # profile -- the class string is exactly the field known to be unreliable.
        self.assertEqual(stats.completeness_profile(4, "Actinopteri"), "core")
        self.assertEqual(stats.completeness_profile(2, "Anthozoa"), "vertebrate")

    def test_class_is_the_fallback_when_no_code_is_given(self):
        self.assertEqual(stats.completeness_profile(None, "Anthozoa"), "core")
        self.assertEqual(stats.completeness_profile(None, "Actinopteri"), "vertebrate")
        self.assertEqual(stats.completeness_profile(None, ""), "vertebrate")


class ProfileAppliedTests(ProcessGffTests):
    """The selected profile actually changes the verdict."""

    def _run_code(self, genes, genetic_code):
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text(gff_for(genes))
        return stats.process_gff(str(p), p.stem, "", 2, genetic_code)

    def test_echinoderm_core_passes_without_the_vertebrate_trnas(self):
        # 13 PCGs + 2 rRNAs, no tRNAs: complete on the core profile.
        core = [g for g in REF if not g.startswith("T") or g.startswith("RNR")]
        s = self._run_code(core, 9)
        self.assertEqual(s["completeness_profile"], "core")
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_correct"], "NA")

    def test_same_annotation_fails_under_the_vertebrate_profile(self):
        core = [g for g in REF if not g.startswith("T") or g.startswith("RNR")]
        s = self._run_code(core, 2)
        self.assertEqual(s["completeness_profile"], "vertebrate")
        self.assertEqual(s["passed"], "no")

    def test_vertebrate_verdict_is_unchanged_by_the_new_argument(self):
        s = self._run_code(REF, 2)
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["missing_genes"], "no")
        self.assertEqual(s["order_correct"], "yes")


if __name__ == "__main__":
    unittest.main()
