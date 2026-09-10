"""Dedup tests for bin/build_invert_reference_db.py.

Both stages exist only because the RefSeq restriction is now per group. While
`refseq_only` is on, a group holds one curated record per genome and these are
no-ops. Lifting it for ctenophora is what makes them load-bearing:

  * NCBI returns a RefSeq record alongside the INSDC submission it was derived
    from (NC_038065 / MG655622, Beroe forskalii) -- near-identical genomes that
    would occupy two of the five slots in a top-n seed panel.
  * It returns NINE Vallicula multiformis isolates. All nine share a family, so
    neither a top-n panel nor select_fallback_seed's family balancing can dilute
    them; only a per-organism cap at build time can.

A rebuild must also be reproducible, so both stages are order-independent.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    HAVE_BIO = True
except ImportError:
    HAVE_BIO = False

if HAVE_BIO:
    import build_invert_reference_db as builder

TWIN_COMMENT = ("PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI "
                "review. The reference sequence is identical to MG655622.1.")


def entry(accession, organism, length, n_cds=11, comment=""):
    rec = SeqRecord(Seq("A" * length), id=accession)
    rec.annotations = {"organism": organism, "comment": comment}
    return (rec, builder.source_of(rec), n_cds, 2)


@unittest.skipUnless(HAVE_BIO, "Biopython not installed")
class RecordSource(unittest.TestCase):
    def test_refseq_is_recognised_by_the_underscore(self):
        self.assertTrue(builder.is_refseq(entry("NC_038065.1", "x", 10)[0]))
        self.assertFalse(builder.is_refseq(entry("MG655622.1", "x", 10)[0]))

    def test_manifest_source_follows_the_record_not_the_search(self):
        """A mixed group must not label every row 'refseq'."""
        self.assertEqual(entry("NC_038065.1", "x", 10)[1], "refseq")
        self.assertEqual(entry("PP327218.1", "x", 10)[1], "genbank")


@unittest.skipUnless(HAVE_BIO, "Biopython not installed")
class DropInsdcTwins(unittest.TestCase):
    def test_drops_the_insdc_copy_and_keeps_the_refseq_one(self):
        kept, n = builder.drop_insdc_twins([
            entry("NC_038065.1", "Beroe forskalii", 13338, comment=TWIN_COMMENT),
            entry("MG655622.1", "Beroe forskalii", 13338),
            entry("PP327218.1", "Tjalfiella sp.", 11397),
        ])
        self.assertEqual(n, 1)
        self.assertEqual([e[0].id for e in kept], ["NC_038065.1", "PP327218.1"])

    def test_keeps_a_different_isolate_of_the_same_species(self):
        """MG655624 is not the twin of NC_038065; only the named accession goes."""
        kept, n = builder.drop_insdc_twins([
            entry("NC_038065.1", "Beroe forskalii", 13338, comment=TWIN_COMMENT),
            entry("MG655624.1", "Beroe forskalii", 13357),
        ])
        self.assertEqual(n, 0)
        self.assertEqual(len(kept), 2)

    def test_is_a_no_op_without_refseq_records(self):
        entries = [entry("PP327218.1", "Tjalfiella sp.", 11397)]
        kept, n = builder.drop_insdc_twins(entries)
        self.assertEqual((kept, n), (entries, 0))


@unittest.skipUnless(HAVE_BIO, "Biopython not installed")
class CapPerOrganism(unittest.TestCase):
    def nine_vallicula(self):
        return [entry(f"PX92269{i}.1", "Vallicula multiformis", 9961 + i) for i in range(9)]

    def test_caps_one_species_without_touching_the_others(self):
        kept, n = builder.cap_per_organism(
            self.nine_vallicula() + [entry("PP327218.1", "Tjalfiella sp.", 11397)], 2)
        self.assertEqual(n, 7)
        self.assertEqual(len(kept), 3)
        self.assertIn("PP327218.1", [e[0].id for e in kept])

    def test_prefers_refseq_then_length(self):
        kept, _n = builder.cap_per_organism([
            entry("MN544300.1", "Hormiphora californensis", 12555),
            entry("NC_045864.1", "Hormiphora californensis", 12564),
            entry("MN544301.1", "Hormiphora californensis", 12000),
        ], 1)
        self.assertEqual([e[0].id for e in kept], ["NC_045864.1"])

    def test_is_order_independent(self):
        entries = self.nine_vallicula()
        self.assertEqual([e[0].id for e in builder.cap_per_organism(entries, 2)[0]],
                         [e[0].id for e in builder.cap_per_organism(list(reversed(entries)), 2)[0]])

    def test_zero_disables_the_cap(self):
        entries = self.nine_vallicula()
        self.assertEqual(builder.cap_per_organism(entries, 0), (entries, 0))


@unittest.skipUnless(HAVE_BIO, "Biopython not installed")
class GroupSpecDefaults(unittest.TestCase):
    def test_only_ctenophora_lifts_the_refseq_restriction(self):
        """Guards the blast radius: lifting it elsewhere re-picks references for
        samples that are already submitted."""
        lifted = {g for g, s in builder.GROUPS.items() if not s.refseq_only}
        self.assertEqual(lifted, {"ctenophora"})

    def test_the_search_term_follows_the_flag(self):
        on = builder.SEARCH_TEMPLATE.format(organism="txid1", refseq=builder.REFSEQ_TERM)
        off = builder.SEARCH_TEMPLATE.format(organism="txid1", refseq="")
        self.assertIn("refseq[filter]", on)
        self.assertNotIn("refseq[filter]", off)
        self.assertIn('mitochondrion[filter] AND "complete genome"[Title]', off)


if __name__ == "__main__":
    unittest.main()
