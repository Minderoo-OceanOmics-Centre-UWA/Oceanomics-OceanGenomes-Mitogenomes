"""The evidence view reports; it never validates.

Its whole reason for existing is that 68 of 68 assemblies in NOVA_260909_LA were
given the identical reason -- "not validated (no Found_in_blast_YN=Yes)" -- for
three quite different situations needing three different actions. The one thing
it must never do is become a second way to pass.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))

from taxonomic_evidence import classify, read_blast

SCRIPT = BIN / "taxonomic_evidence.py"


def blast_row(species, pident, aln, qlen):
    f = [""] * 9
    f[4], f[6], f[7], f[8] = species, str(pident), str(aln), str(qlen)
    return "\t".join(f)


class NeverValidatesTests(unittest.TestCase):
    def test_no_write_path_to_any_authoritative_store(self):
        """Structural, not behavioural: the script must not import psycopg2, and
        must not mention lca_validation or Found_in_blast_YN as something it sets."""
        text = SCRIPT.read_text()
        self.assertNotIn("import psycopg2", text)
        self.assertNotIn("upsert_lca_validation", text)
        self.assertNotIn("INSERT INTO", text)

    def test_validated_flag_is_read_not_decided(self):
        """nominal_id_validated mirrors the authoritative result; strong evidence
        alone must never flip it."""
        hits = [("Chaceon granulatus", 99.9, 1.0)]
        _r, _t, outcome, _q, _a = classify("Chaceon sp.", {"genus": "Chaceon"},
                                           hits, [], validated=False)
        self.assertNotEqual(outcome, "nominal_id_validated")


class OutcomeTests(unittest.TestCase):
    def test_above_species_label_is_not_a_failure(self):
        """61 of 68 assemblies land here. 'Never identified to species' is a
        different instruction from 'the assembly contradicts the label'."""
        _r, _t, outcome, _q, _a = classify(
            "Tetractinellida (Order)", {"order": "Tetractinellida"}, [], [], False)
        self.assertEqual(outcome, "no_species_level_claim")

    def test_rank_is_tested_before_hit_availability(self):
        """A label above species can never validate no matter what the databases
        hold, so that is the more useful answer even when there are no hits."""
        hits = [("Something else", 99.0, 1.0)]
        _r, _t, outcome, _q, _a = classify(
            "Hexactinellida (Class)", {"class": "Hexactinellida"}, hits, [], False)
        self.assertEqual(outcome, "no_species_level_claim")

    def test_discordant_high_identity_is_called_out(self):
        """A species claim contradicted at high identity is a specimen labelling
        problem, not an assembly problem, and should read differently."""
        hits = [("Psenes cyanophrys", 100.0, 1.0)]
        _r, _t, outcome, _q, _a = classify("Ariomma indicum", {}, hits, [], False)
        self.assertEqual(outcome, "high_identity_evidence_discordant")

    def test_unresolved_label_is_its_own_outcome(self):
        _r, _t, outcome, _q, _a = classify("Larval fish", {}, [], [], False)
        self.assertEqual(outcome, "taxonomy_unresolved")

    def test_validated_sample_reports_validated(self):
        _r, _t, outcome, _q, _a = classify("Chaceon sp.", {}, [], [], True)
        self.assertEqual(outcome, "nominal_id_validated")


class SubThresholdEvidenceTests(unittest.TestCase):
    def test_sub_threshold_hits_are_kept_with_their_span(self):
        """OG3060's three ~98.5% octocoral hits span 20% of the query, so the
        validator correctly rejects them -- and the run currently discards them
        entirely, which is what leaves a curator with nothing to act on."""
        with tempfile.TemporaryDirectory() as tmp:
            blast = Path(tmp) / "b.tsv"
            blast.write_text("\n".join([
                blast_row("Chrysogorgia tricaulis", 98.588, 354, 1788),
                blast_row("Ideogorgia capensis", 98.305, 354, 1788),
            ]) + "\n")
            out = Path(tmp) / "out.tsv"
            subprocess.run(
                [sys.executable, str(SCRIPT), "--assembly", "A", "--og-id", "OG3060",
                 "--nominal", "Scleralcyonacea sp.", "--blast", str(blast),
                 "--order", "Scleralcyonacea", "--out", str(out)],
                check=True, capture_output=True,
                env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"})
            body = out.read_text()
        self.assertIn("Chrysogorgia tricaulis", body)
        self.assertIn("nominal_id_validated", body.splitlines()[0])
        row = dict(zip(body.splitlines()[0].split("\t"), body.splitlines()[1].split("\t")))
        self.assertEqual(row["nominal_id_validated"], "false")
        self.assertEqual(row["best_hit_qualifying"], "false")
        self.assertNotEqual(row["sub_threshold_hits"], "none")

    def test_query_cover_is_computed_from_alignment_and_query_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "b.tsv"
            p.write_text(blast_row("X", 98.5, 354, 1788) + "\n")
            hits = read_blast(p)
        self.assertAlmostEqual(hits[0][2], 354 / 1788, places=4)


if __name__ == "__main__":
    unittest.main()
