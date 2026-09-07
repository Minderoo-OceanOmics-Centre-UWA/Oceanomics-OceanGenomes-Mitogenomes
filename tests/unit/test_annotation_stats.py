"""Unit tests for bin/annotation_stats.py -- the completeness gate.

Covers process_gff(): missing_genes / order_correct / passed and the bounded
tRNA-tolerance branch (trna_advisory), plus the reduced-tRNA classification.
Cnidaria (corals, anemones, jellyfish) and Porifera (sponges) both get the
relaxed completeness bar -- see InvertTaxonGroups in lib/ for why these two
groups are handled together while the rest of the invertebrate batch is not.
Pure stdlib -- the Biopython import lives inside process_protein_lengths(),
which these tests do not call -- so they always run.
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


def gff_for(genes, region_len=26500, spacing=0, gaps=None):
    """Minimal EMMA-style GFF: one `gene` line + one matching feature line per
    name, at ascending coordinates in the given order.

    `spacing` is the intergenic gap in bp between every consecutive pair, and
    `gaps` overrides it for individual pairs, keyed by the index of the gene on
    the LEFT of the pair.

    Spacing is explicit and defaults to abutting genes ON PURPOSE. It used to be
    implicit: genes were laid out as pos..pos+50 stepping pos += 100, so every
    consecutive pair sat exactly 49 bp apart and every test in this file cleared
    the 50 bp annotation_gaps default by a single base. That was luck, not design
    -- a later change to the threshold would have silently made thirty-six phantom
    gaps appear across the whole suite.
    """
    gaps = gaps or {}
    lines = ["##gff-version 3", f"##sequence-region chr 1 {region_len}"]
    pos = 100
    for idx, g in enumerate(genes):
        end = pos + 50
        gid = f"gene-{g}"
        lines.append(f"chr\tEmma\tgene\t{pos}\t{end}\t.\t+\t.\tID={gid};Name=MT-{g}")
        lines.append(
            f"chr\tEmma\t{_feature_type(g)}\t{pos}\t{end}\t.\t+\t.\t"
            f"ID=feat-{g};Parent={gid};Name=MT-{g}")
        pos = end + gaps.get(idx, spacing) + 1
    return "\n".join(lines) + "\n"


def order_with_swap(first, second):
    """REF_GENES with the adjacent pair (first, second) transposed."""
    genes = list(REF)
    i = genes.index(first)
    assert genes[i + 1] == second, f"{first},{second} are not adjacent in REF_GENES"
    genes[i:i + 2] = [second, first]
    return genes


class ProcessGffTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _run(self, genes, trna_tolerance=2, class_name="", taxon=None,
             genetic_code=None, spacing=0, gaps=None, gap_threshold=50):
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text(gff_for(genes, spacing=spacing, gaps=gaps))
        return stats.process_gff(str(p), p.stem, class_name, trna_tolerance,
                                 genetic_code, taxon=taxon,
                                 gap_threshold=gap_threshold)

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



class ReducedTrnaExpectationTests(unittest.TestCase):
    def test_cnidaria_classes_are_reduced_trna(self):
        for class_name in ["Anthozoa", "Hydrozoa", "Scyphozoa", "Cnidaria"]:
            self.assertTrue(stats.has_reduced_trna_expectation(class_name))

    def test_porifera_classes_are_reduced_trna(self):
        for class_name in ["Demospongiae", "Calcarea", "Hexactinellida", "Porifera"]:
            self.assertTrue(stats.has_reduced_trna_expectation(class_name))

    def test_case_and_whitespace_insensitive(self):
        self.assertTrue(stats.has_reduced_trna_expectation("  porifera  "))
        self.assertTrue(stats.has_reduced_trna_expectation("ANTHOZOA"))

    def test_other_invert_phyla_are_not_reduced_trna(self):
        for class_name in ["Gastropoda", "Bivalvia", "Malacostraca", "Pycnogonida",
                            "Asteroidea", "Echinoidea", "Actinopteri", None, ""]:
            self.assertFalse(stats.has_reduced_trna_expectation(class_name))


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


class OrderVariantTests(unittest.TestCase):
    """The curated per-taxon gene-order table.

    The gate used to test one hard-coded vertebrate order, so a real published
    lineage-level rearrangement failed exactly like a scrambled assembly. It now
    consults bin/mito_gene_order.py ORDER_VARIANTS, keyed at the rank the evidence
    supports, and records which rule fired.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _run(self, genes, taxon=None, gaps=None, gap_threshold=50):
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text(gff_for(genes, gaps=gaps))
        return stats.process_gff(str(p), p.stem, genetic_code=2, taxon=taxon,
                                 gap_threshold=gap_threshold)

    # -- the variant is ACCEPTED, not REQUIRED ------------------------------

    def test_canonical_order_in_a_variant_taxon_still_passes(self):
        # A rule ADDS an accepted order; it does not replace the canonical one.
        # A canonical member of a rule-carrying taxon must keep passing and must
        # NOT acquire an order_variant value merely because its genus has a rule.
        s = self._run(REF, taxon={"genus": "Chlorurus"})
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_correct"], "yes")
        self.assertEqual(s["order_variant"], "no")

    def test_imq_order_with_the_keyed_genus_passes_as_a_variant(self):
        imq = order_with_swap("TQ", "TM")
        for genus in ("Chlorurus", "Hipposcarus"):
            s = self._run(imq, taxon={"genus": genus})
            self.assertEqual(s["passed"], "yes", genus)
            self.assertEqual(s["order_correct"], "variant", genus)
            self.assertEqual(s["order_variant"], "scarine_imq", genus)

    def test_matching_is_case_and_whitespace_insensitive(self):
        imq = order_with_swap("TQ", "TM")
        s = self._run(imq, taxon={"genus": "  chlorurus  "})
        self.assertEqual(s["order_variant"], "scarine_imq")

    def test_diploprion_ds_variant(self):
        ds = order_with_swap("TS2", "TD")
        s = self._run(ds, taxon={"genus": "Diploprion"})
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_variant"], "diploprion_ds")

    # -- the rule must not leak --------------------------------------------

    def test_the_same_order_in_an_unkeyed_genus_still_fails(self):
        # Epibulus is Labridae and canonically IQM, while the parrotfishes (also
        # Labridae under current classification) are IMQ -- which is exactly why
        # these rules are genus-keyed and a family key for Labridae would be wrong.
        imq = order_with_swap("TQ", "TM")
        s = self._run(imq, taxon={"genus": "Epibulus"})
        self.assertEqual(s["passed"], "no")
        self.assertEqual(s["order_correct"], "no")
        self.assertEqual(s["order_variant"], "no")

    def test_the_rule_does_not_leak_across_ranks(self):
        imq = order_with_swap("TQ", "TM")
        for taxon in ({"order": "Perciformes"}, {"family": "Scaridae"},
                      {"class": "Actinopteri"}):
            s = self._run(imq, taxon=taxon)
            self.assertEqual(s["passed"], "no", taxon)

    def test_no_taxonomy_at_all_fails(self):
        # '', None, [] and the literal 'unknown' all occur in real samplesheets.
        imq = order_with_swap("TQ", "TM")
        for taxon in (None, {}, {"genus": ""}, {"genus": None}, {"genus": []},
                      {"genus": "unknown"}, {"genus": "NA"},
                      {"genus": "", "family": "unknown", "order": None, "class": []}):
            s = self._run(imq, taxon=taxon)
            self.assertEqual(s["passed"], "no", taxon)
            self.assertEqual(s["order_variant"], "no", taxon)

    def test_a_scramble_in_a_keyed_genus_still_fails(self):
        genes = list(REF)
        genes[5], genes[25] = genes[25], genes[5]
        s = self._run(genes, taxon={"genus": "Chlorurus"})
        self.assertEqual(s["passed"], "no")

    # -- rank precedence ---------------------------------------------------

    def test_most_specific_rank_wins_outright_and_ranks_do_not_compose(self):
        # Injected rows, so the behaviour is pinned even though no live rule pair
        # currently exercises it -- and an ORDER-ranked key is tested at all, since
        # RANK_PRECEDENCE names that rank and it is where the next rule may land.
        mgo = sys.modules.get("mito_gene_order") or __import__("mito_gene_order")
        family_key = ("family", "Testidae")
        order_key = ("order", "Testiformes")
        mgo.ORDER_VARIANTS[family_key] = [
            ("family_rule", ("TQ", "TM"), ("TM", "TQ"), "test")]
        mgo.ORDER_VARIANTS[order_key] = [
            ("order_rule", ("TS2", "TD"), ("TD", "TS2"), "test")]
        try:
            # Family beats order.
            s = self._run(order_with_swap("TQ", "TM"),
                          taxon={"family": "Testidae", "order": "Testiformes"})
            self.assertEqual(s["order_variant"], "family_rule")
            # And the order rule does NOT also apply: ranks do not compose.
            s = self._run(order_with_swap("TS2", "TD"),
                          taxon={"family": "Testidae", "order": "Testiformes"})
            self.assertEqual(s["passed"], "no")
            # With no family, the order-ranked key does fire.
            s = self._run(order_with_swap("TS2", "TD"), taxon={"order": "Testiformes"})
            self.assertEqual(s["order_variant"], "order_rule")
        finally:
            del mgo.ORDER_VARIANTS[family_key]
            del mgo.ORDER_VARIANTS[order_key]

    # -- the two published clade rules --------------------------------------

    def _anguilliform(self):
        """ND6+trnE moved from upstream of CYTB to between trnT and trnP.

        Published as trnT -> [control region] -> ND6 -> trnE -> trnP; since the
        control region is not a gene and REF_GENES starts at trnF, that renders as
        this order.
        """
        genes = [g for g in REF if g not in ("ND6", "TE")]
        i = genes.index("TT") + 1
        return genes[:i] + ["ND6", "TE"] + genes[i:]

    def _macrourid(self):
        """trnE alone moved to after trnP: the published trnT-trnP-trnE cluster."""
        genes = [g for g in REF if g != "TE"]
        return genes + ["TE"]

    def test_anguilliform_order_passes_for_each_keyed_family(self):
        for family in ("Congridae", "Nettastomatidae", "Colocongridae",
                       "Muraenesocidae"):
            s = self._run(self._anguilliform(), taxon={"family": family})
            self.assertEqual(s["passed"], "yes", family)
            self.assertEqual(s["order_variant"], "anguilliform_nd6_te", family)

    def test_anguilliform_order_fails_in_the_families_that_are_canonical(self):
        # THE reason this rule is family-keyed and not order-keyed. Anguillidae,
        # Synaphobranchidae, Muraenidae and Serrivomeridae are published as
        # retaining the typical vertebrate order; an order-rank key would grant
        # licence across all four. If either of these passes, the rule has been
        # widened back to the order rank.
        for family in ("Synaphobranchidae", "Nemichthyidae", "Anguillidae",
                       "Muraenidae", "Serrivomeridae"):
            s = self._run(self._anguilliform(),
                          taxon={"family": family, "order": "Anguilliformes"})
            self.assertEqual(s["passed"], "no", family)

    def test_a_canonical_anguilliform_keeps_passing_with_no_variant(self):
        for family in ("Congridae", "Synaphobranchidae"):
            s = self._run(REF, taxon={"family": family})
            self.assertEqual(s["passed"], "yes", family)
            self.assertEqual(s["order_variant"], "no", family)

    def test_the_blachea_genus_row_reaches_what_the_family_row_cannot(self):
        # Blachea is Colocongridae but resolves to a blank family in the
        # samplesheet, so the family row cannot reach it.
        s = self._run(self._anguilliform(),
                      taxon={"genus": "Blachea", "family": "", "order": ""})
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_variant"], "anguilliform_nd6_te")

    def test_the_blachea_row_covers_only_the_genus_it_names(self):
        # The paired negative: another congrid genus with the same blank taxonomy
        # gets no rule, which is a deliberate choice rather than a surprise.
        s = self._run(self._anguilliform(),
                      taxon={"genus": "Ariosoma", "family": "", "order": ""})
        self.assertEqual(s["passed"], "no")

    def test_macrourid_order_passes_for_macrouridae(self):
        s = self._run(self._macrourid(), taxon={"family": "Macrouridae"})
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_variant"], "macrourid_te")

    def test_macrourid_rule_does_not_leak_to_other_gadiforms(self):
        s = self._run(self._macrourid(), taxon={"family": "Gadidae"})
        self.assertEqual(s["passed"], "no")

    def test_a_canonical_macrourid_keeps_passing_without_a_variant(self):
        # The Bathygadus case. A canonical member of a rule-carrying family must
        # not acquire an order_variant merely because its family has a rule --
        # which matters here because Macrouridae carries several DIFFERENT
        # rearrangement patterns across its subfamilies.
        s = self._run(REF, taxon={"family": "Macrouridae"})
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_correct"], "yes")
        self.assertEqual(s["order_variant"], "no")

    def test_the_two_rules_are_distinct(self):
        # macrourid_te moves trnE alone; anguilliform_nd6_te moves ND6 with it.
        # Neither taxon may accept the other's order.
        self.assertEqual(
            self._run(self._macrourid(), taxon={"family": "Congridae"})["passed"], "no")
        self.assertEqual(
            self._run(self._anguilliform(), taxon={"family": "Macrouridae"})["passed"], "no")

    def test_a_malformed_table_entry_raises(self):
        mgo = sys.modules.get("mito_gene_order") or __import__("mito_gene_order")
        bad_key = ("genus", "Malformed")
        for rule, why in [
            (("nc", ("TQ", "CO1"), ("CO1", "TQ"), "x"), "non-contiguous block"),
            (("np", ("TQ", "TM"), ("TQ", "TF"), "x"), "not a permutation"),
        ]:
            mgo.ORDER_VARIANTS[bad_key] = [rule]
            try:
                with self.assertRaises(ValueError, msg=why):
                    mgo.ref_order_for({"genus": "Malformed"})
            finally:
                del mgo.ORDER_VARIANTS[bad_key]


class AnnotationGapTests(unittest.TestCase):
    """Interior intergenic spans. Purely advisory: they never affect `passed`.

    This is what separates a real transposition (a tRNA moved) from a missed call
    (a tRNA-shaped hole left at both the origin and the destination of the
    apparent move) -- a distinction otherwise only visible by reading the GFF.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _run(self, genes, gaps=None, gap_threshold=50, genetic_code=2, taxon=None,
             core_gap_threshold=stats.DEFAULT_CORE_GAP_THRESHOLD):
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text(gff_for(genes, region_len=80000, gaps=gaps))
        return stats.process_gff(str(p), p.stem, genetic_code=genetic_code,
                                 taxon=taxon, gap_threshold=gap_threshold,
                                 core_gap_threshold=core_gap_threshold)

    def test_a_large_interior_gap_is_reported_and_passed_is_unchanged(self):
        i = REF.index("ND4")
        s = self._run(REF, gaps={i: 213})
        self.assertIn("ND4:", s["annotation_gaps"])
        self.assertIn("(213)", s["annotation_gaps"])
        self.assertEqual(s["passed"], "yes")
        self.assertEqual(s["order_correct"], "yes")

    def test_abutting_genes_report_nothing(self):
        self.assertEqual(self._run(REF)["annotation_gaps"], "no")

    def test_the_control_region_is_exempt(self):
        # TP -> TF is the vertebrate control region and is legitimately ~1 kb.
        genes = ["TP", "TF"] + [g for g in REF if g not in ("TP", "TF")]
        self.assertEqual(self._run(genes, gaps={0: 1100})["annotation_gaps"], "no")

    def test_threshold_boundary(self):
        # 62 bp is the missed-tRNA signal the default has to reach; 35 bp is an
        # OriL-sized span that must stay silent; 51 bp is also OriL-sized and
        # DOES fire, because the real OriL range is 30-51 bp and straddles the
        # boundary. That is expected, and is the price of reaching 62.
        i = REF.index("ND4")
        self.assertEqual(self._run(REF, gaps={i: 35})["annotation_gaps"], "no")
        self.assertIn("(51)", self._run(REF, gaps={i: 51})["annotation_gaps"])
        self.assertIn("(62)", self._run(REF, gaps={i: 62})["annotation_gaps"])

    def test_the_threshold_is_configurable(self):
        i = REF.index("ND4")
        self.assertEqual(self._run(REF, gaps={i: 62}, gap_threshold=100)["annotation_gaps"], "no")

    def test_gaps_are_reported_on_the_core_profile_at_its_own_threshold(self):
        # The invert-facing half. A core-profile assembly is judged on gene
        # presence alone with order reported NA, so without this a coral with all
        # its core genes and 30 kb of unannotated sequence between two of them
        # passes with nothing recorded at all.
        i = REF.index("ND4")
        for code in (4, 5, 9):
            s = self._run(REF, gaps={i: 30832}, genetic_code=code)
            self.assertEqual(s["completeness_profile"], "core", code)
            self.assertEqual(s["order_correct"], "NA", code)
            self.assertIn("(30832)", s["annotation_gaps"], code)
            self.assertEqual(s["passed"], "yes", code)

    def test_the_core_profile_uses_its_own_much_higher_threshold(self):
        # The two profiles look for different things and one threshold cannot
        # serve both. Measured on batch-20's 38 coral assemblies, the vertebrate
        # 50 bp default reports 7.3 gaps per assembly and flags every one of them.
        # A 213 bp gap is a real signal on a vertebrate and unremarkable on a coral.
        i = REF.index("ND4")
        self.assertIn("(213)", self._run(REF, gaps={i: 213}, genetic_code=2)["annotation_gaps"])
        self.assertEqual(self._run(REF, gaps={i: 213}, genetic_code=4)["annotation_gaps"], "no")

    def test_the_conserved_coral_intergenic_region_stays_silent(self):
        # Batch-20's corals share a 1725 bp CO3->CO2 span, byte-identical across
        # 20 of 38 assemblies -- plainly biology, not a defect. At any threshold
        # low enough to report it, that one span alone flags half the batch.
        i = REF.index("CO3")
        self.assertEqual(self._run(REF, gaps={i: 1725}, genetic_code=4)["annotation_gaps"], "no")

    def test_the_core_threshold_is_configurable(self):
        i = REF.index("ND4")
        s = self._run(REF, gaps={i: 213}, genetic_code=4, core_gap_threshold=100)
        self.assertIn("(213)", s["annotation_gaps"])


class OrderDeviationTests(unittest.TestCase):
    """Triage for the held pile. Advisory: never affects `passed`."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _run(self, genes, taxon=None):
        p = Path(self.tmp) / "OG1.ilmn.240101.getorg1770.emma102.gff"
        p.write_text(gff_for(genes))
        return stats.process_gff(str(p), p.stem, genetic_code=2, taxon=taxon)

    def test_correct_order_reports_nothing(self):
        self.assertEqual(self._run(REF)["order_deviation"], "no")

    def test_a_matched_variant_reports_nothing(self):
        s = self._run(order_with_swap("TQ", "TM"), taxon={"genus": "Chlorurus"})
        self.assertEqual(s["order_deviation"], "no")

    def test_one_displaced_gene_scores_low(self):
        s = self._run(order_with_swap("TQ", "TM"))
        self.assertEqual(s["order_deviation"], "1")
        self.assertEqual(s["passed"], "no")

    def test_a_reversed_segment_scores_much_higher(self):
        s = self._run(REF[:18] + REF[18:][::-1])
        self.assertGreater(int(s["order_deviation"]), 10)
        self.assertEqual(s["passed"], "no")

    def test_deviation_is_measured_against_the_nearest_accepted_order(self):
        # A taxon with a rule is scored against whichever ordering it is closer
        # to, rather than being penalised for the rule existing.
        genes = order_with_swap("TQ", "TM")
        with_rule = self._run(genes, taxon={"genus": "Chlorurus"})
        self.assertEqual(with_rule["order_deviation"], "no")  # it matched outright
        # One further swap on top of the variant: still near the variant order.
        genes2 = list(genes)
        i = genes2.index("TS2")
        genes2[i:i + 2] = [genes2[i + 1], genes2[i]]
        s = self._run(genes2, taxon={"genus": "Chlorurus"})
        self.assertEqual(s["passed"], "no")
        self.assertEqual(s["order_deviation"], "1")
