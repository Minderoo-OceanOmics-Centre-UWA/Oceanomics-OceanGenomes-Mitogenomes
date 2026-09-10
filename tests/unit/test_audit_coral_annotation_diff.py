"""Logic tests for bin/audit_coral_annotation_diff.py.

The audit re-annotates published corals against a rebuilt reference database, so its
answer decides whether a database rebuild is safe for the one consumer that copies
CONTENT out of the reference. Two things have to be right or the answer is worse than
useless:

  * The frame. The cox1-rotated genome the MITOS BED refers to is not published, so it is
    reconstructed. annotation/<prefix>.fa looks like it but is the published RE-ORIGINED
    genome, a different rotation. Annotating against the wrong one produces coordinates
    that are confidently wrong rather than obviously broken, so a sample whose frame
    cannot be verified must fail, never proceed. Two earlier designs were wrong against
    real data: "cox1 sits at offset 0" rejects rotate_to_cox1.py's unrotated fail-safe,
    and "the BED cox1 matches the published CO1" rejects any coral whose cox1 the fixer
    repaired. The check therefore compares genes the fixer never touches. OG2361 is the
    real sample that caught both.
  * The attribution. The published annotation was produced by older code on this branch,
    so a difference between it and a fresh run is not evidence about the reference. The
    third arm (today's code, OLD reference) is what separates code drift from a
    reference-driven change, and a difference that shows up there must never be blamed
    on the reference.

Everything under test is pure, so it runs without biopython, BLAST or a database.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

import audit_coral_annotation_diff as audit  # noqa: E402

# A scleractinian BED as ROTATE_ORIGIN leaves it: cox1 first, starting at offset 0.
ROTATED_BED = [
    ("OG2358", 0, 1602, "cox1", "+"),
    ("OG2358", 1796, 1867, "trnM(cat)", "+"),
    ("OG2358", 2076, 4025, "rrnL", "+"),
    ("OG2358", 4360, 5070, "nad5_0", "+"),
]
# OG2361: rotate_to_cox1.py found no confident cox1 hit and wrote the assembly through
# UNROTATED, so MITOS annotated that frame. cox1 is a fragment beside a group-I intron
# and nothing starts at 0. This reconstruction is CORRECT and must be accepted.
UNROTATED_BED = [
    ("OG2361", 30, 729, "OH", "+"),
    ("OG2361", 719, 1051, "gpI", "+"),
    ("OG2361", 1750, 2623, "cox1", "+"),
]

# A genome whose cox1 interval holds something else entirely: a reconstruction that
# landed in the wrong frame.
WRONG_FRAME = "T" * 30000


def metrics(sha, pcgs=13, nad5=605, co1=533, stops=0, failed=False):
    return {"sha": sha, "pcgs": pcgs, "nad5_aa": nad5, "co1_aa": co1,
            "internal_stops": stops, "failed": failed}


UNTOUCHED = {"CO2": "ATG" + "GCT" * 100, "CO3": "ATG" + "TTC" * 90,
             "CYTB": "ATG" + "AAG" * 110, "ND4": "ATG" + "CCA" * 120}


def genome_from(bed, genes, total=20000):
    """A synthetic genome carrying each gene's sequence at its BED coordinates."""
    seq = list("A" * total)
    for _chrom, start, end, name, _strand in bed:
        gene = audit.emma_gene(name)
        if gene in genes:
            piece = genes[gene]
            seq[start:start + len(piece)] = list(piece)
    return "".join(seq)


# A BED holding four genes the coral fixer never rewrites, plus the three it does.
MIXED_BED = [
    ("OG", 100, 100 + len(UNTOUCHED["CO2"]), "cox2", "+"),
    ("OG", 2000, 2000 + len(UNTOUCHED["CO3"]), "cox3", "+"),
    ("OG", 4000, 4000 + len(UNTOUCHED["CYTB"]), "cob", "+"),
    ("OG", 6000, 6000 + len(UNTOUCHED["ND4"]), "nad4", "+"),
    ("OG", 8000, 8873, "cox1", "+"),      # fixer-repaired, must be ignored
    ("OG", 10000, 11000, "nad5_0", "+"),  # fixer-repaired, must be ignored
]


class GeneMapping(unittest.TestCase):
    def test_mitos_names_map_to_emma_names(self):
        self.assertEqual(audit.emma_gene("cox1"), "CO1")
        self.assertEqual(audit.emma_gene("nad4l"), "ND4L")
        self.assertEqual(audit.emma_gene("nad5_0"), "ND5")
        self.assertEqual(audit.emma_gene("rrnL"), "RNR2")
        self.assertEqual(audit.emma_gene("trnM(cat)"), "")
        self.assertEqual(audit.emma_gene(""), "")

    def test_the_copied_map_has_not_drifted_from_mitos_to_emma(self):
        """The map is copied rather than imported (mitos_to_emma needs biopython at
        import). This is what stops the copy going stale."""
        try:
            import mitos_to_emma
        except ImportError:
            self.skipTest("biopython not installed")
        self.assertEqual(audit.PCG_MAP, mitos_to_emma.PCG_MAP)
        self.assertEqual(audit.RRNA_MAP, mitos_to_emma.RRNA_MAP)


class FrameCheck(unittest.TestCase):
    def test_a_matching_frame_passes_on_untouched_genes(self):
        genome = genome_from(MIXED_BED, UNTOUCHED)
        ok, reason = audit.check_frame(MIXED_BED, genome, UNTOUCHED)
        self.assertTrue(ok, reason)
        self.assertIn("4/4", reason)

    def test_a_repaired_cox1_does_not_break_the_frame_check(self):
        """OG2361: raw-BED cox1 is 873 bp while the published CO1 is 1572 bp because the
        fixer repaired it. Anchoring on cox1 rejected a correct reconstruction."""
        published = dict(UNTOUCHED, CO1="ATG" + "GGG" * 523, ND5="ATG" + "TTT" * 200)
        genome = genome_from(MIXED_BED, UNTOUCHED)
        ok, reason = audit.check_frame(MIXED_BED, genome, published)
        self.assertTrue(ok, reason)

    def test_a_shifted_frame_is_rejected(self):
        genome = genome_from(MIXED_BED, UNTOUCHED)
        ok, reason = audit.check_frame(MIXED_BED, "TTTT" + genome, UNTOUCHED)
        self.assertFalse(ok)
        self.assertIn("does not match", reason)

    def test_a_frame_with_too_few_agreements_is_rejected(self):
        """One accidental agreement is not evidence. Requiring several independent genes
        is what makes a content check strong."""
        genome = genome_from(MIXED_BED, {"CO2": UNTOUCHED["CO2"]})
        ok, reason = audit.check_frame(MIXED_BED, genome, UNTOUCHED)
        self.assertFalse(ok)
        self.assertIn("1/4", reason)

    def test_a_genome_shorter_than_the_coordinates_is_rejected(self):
        ok, reason = audit.check_frame(MIXED_BED, "A" * 500, UNTOUCHED)
        self.assertFalse(ok)
        self.assertIn("exceeds genome length", reason)

    def test_an_empty_reconstruction_is_rejected(self):
        ok, _ = audit.check_frame(MIXED_BED, "", UNTOUCHED)
        self.assertFalse(ok)

    def test_only_fixer_touched_genes_available_cannot_anchor_a_frame(self):
        bed = [("OG", 0, 873, "cox1", "+"), ("OG", 1000, 2000, "nad5_0", "+")]
        ok, reason = audit.check_frame(bed, "A" * 5000, {"CO1": "AAA", "ND5": "AAA"})
        self.assertFalse(ok)
        self.assertIn("no fixer-independent gene", reason)

    def test_minus_strand_features_are_reverse_complemented(self):
        gene = UNTOUCHED["CO2"]
        bed = [("OG", 100, 100 + len(gene), "cox2", "-"),
               ("OG", 2000, 2000 + len(UNTOUCHED["CO3"]), "cox3", "+"),
               ("OG", 4000, 4000 + len(UNTOUCHED["CYTB"]), "cob", "+"),
               ("OG", 6000, 6000 + len(UNTOUCHED["ND4"]), "nad4", "+")]
        genome = list(genome_from(bed, UNTOUCHED))
        genome[100:100 + len(gene)] = list(audit.revcomp(gene))
        ok, reason = audit.check_frame(bed, "".join(genome), UNTOUCHED)
        self.assertTrue(ok, reason)
        self.assertIn("4/4", reason)


class Attribution(unittest.TestCase):
    """published -> oldref isolates code drift; oldref -> newref isolates the reference."""

    def test_a_reference_driven_change_is_attributed_to_the_reference(self):
        verdict, attribution = audit.verdict_of(metrics("aaa"), metrics("aaa"),
                                                metrics("bbb"))
        self.assertEqual(attribution, "reference")
        self.assertNotEqual(verdict, "identical")

    def test_code_drift_is_never_blamed_on_the_reference(self):
        """The whole reason the third arm exists. The published output differs from a
        fresh run, but both references produce the same annotation, so the reference
        changed nothing."""
        verdict, attribution = audit.verdict_of(metrics("aaa"), metrics("bbb"),
                                                metrics("bbb"))
        self.assertEqual(attribution, "code_drift")
        self.assertNotEqual(verdict, "identical")

    def test_both_moving_is_reported_as_both(self):
        _verdict, attribution = audit.verdict_of(metrics("aaa"), metrics("bbb"),
                                                 metrics("ccc"))
        self.assertEqual(attribution, "both")

    def test_nothing_moving_is_identical_with_no_attribution(self):
        verdict, attribution = audit.verdict_of(metrics("aaa"), metrics("aaa"),
                                                metrics("aaa"))
        self.assertEqual((verdict, attribution), ("identical", "n/a"))

    def test_an_unchanged_reference_leaves_the_two_reruns_identical(self):
        """The built-in control: 27 of the batch-20 corals keep the same reference, so
        their oldref and newref arms must agree. If they do not, the harness is
        nondeterministic and no other verdict can be trusted."""
        _verdict, attribution = audit.verdict_of(metrics("aaa"), metrics("bbb"),
                                                 metrics("bbb"))
        self.assertNotEqual(attribution, "reference")

    def test_a_failed_arm_is_failed_not_silently_identical(self):
        verdict, attribution = audit.verdict_of(metrics("aaa"),
                                                metrics("aaa", failed=True),
                                                metrics("aaa"))
        self.assertEqual((verdict, attribution), ("failed", "n/a"))

    def test_a_missing_arm_is_failed(self):
        self.assertEqual(audit.verdict_of(metrics("aaa"), None, metrics("aaa")),
                         ("failed", "n/a"))


class QualityDirection(unittest.TestCase):
    def test_recovering_a_pcg_reads_as_better(self):
        verdict, _ = audit.verdict_of(metrics("aaa", pcgs=12), metrics("aaa", pcgs=12),
                                      metrics("bbb", pcgs=13))
        self.assertEqual(verdict, "changed_better")

    def test_losing_a_pcg_reads_as_worse(self):
        verdict, _ = audit.verdict_of(metrics("aaa", pcgs=13), metrics("aaa", pcgs=13),
                                      metrics("bbb", pcgs=12))
        self.assertEqual(verdict, "changed_worse")

    def test_an_introduced_internal_stop_outweighs_a_longer_nad5(self):
        """A broken reading frame is not compensated by length: an internal stop is the
        defect a transferred nad5 join can introduce, and it must dominate."""
        verdict, _ = audit.verdict_of(metrics("aaa", nad5=605, stops=0),
                                      metrics("aaa", nad5=605, stops=0),
                                      metrics("bbb", nad5=640, stops=1))
        self.assertEqual(verdict, "changed_worse")

    def test_a_different_sequence_at_equal_quality_is_neutral(self):
        verdict, _ = audit.verdict_of(metrics("aaa"), metrics("aaa"), metrics("bbb"))
        self.assertEqual(verdict, "changed_neutral")

    def test_internal_stops_ignore_the_terminal_stop(self):
        self.assertEqual(audit.count_internal_stops(["MKV*"]), 0)
        self.assertEqual(audit.count_internal_stops(["MK*V*"]), 1)
        self.assertEqual(audit.count_internal_stops(["MKV"]), 0)
        self.assertEqual(audit.count_internal_stops(["M*K*V*", "MKV*"]), 2)


class Circularity(unittest.TestCase):
    """meta.circular == false is what makes the module pass --linear, and a linear genome
    must not be re-origined as if it wrapped."""

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)

    def write(self, verdict):
        p = self.dir / "check.tsv"
        p.write_text("sample\tgetorg_circular\tfinal_verdict_circular\n"
                     f"OG2358\tTrue\t{verdict}\n")
        return p

    def test_true_is_circular(self):
        self.assertTrue(audit.circular_from_getorg_check(self.write("True")))

    def test_false_is_linear(self):
        self.assertFalse(audit.circular_from_getorg_check(self.write("False")))

    def test_na_defaults_to_circular(self):
        self.assertTrue(audit.circular_from_getorg_check(self.write("NA")))

    def test_a_missing_file_defaults_to_circular(self):
        self.assertTrue(audit.circular_from_getorg_check(self.dir / "absent.tsv"))


class ReferenceParsing(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)

    def test_accession_and_organism_from_a_genbank_header(self):
        gb = self.dir / "ref.gb"
        gb.write_text("LOCUS       NC_040137  17887 bp    DNA     circular INV\n"
                      "ACCESSION   NC_040137\n"
                      "VERSION     NC_040137.1\n"
                      "  ORGANISM  Montipora efflorescens\n"
                      "            Eukaryota; Metazoa; Cnidaria.\n")
        self.assertEqual(audit.reference_accession(gb),
                         ("NC_040137.1", "Montipora efflorescens"))

    def test_a_missing_file_is_blank_not_an_exception(self):
        self.assertEqual(audit.reference_accession(self.dir / "absent.gb"), ("", ""))


if __name__ == "__main__":
    unittest.main()
