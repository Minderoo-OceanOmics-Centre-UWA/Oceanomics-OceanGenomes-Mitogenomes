import configparser
import csv
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "species_validation", ROOT / "bin" / "species_validation.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FakeCursor:
    """Routes fetchone() by the most recent query so one fake DB can back both
    get_species_for_ogid (SELECT ... FROM sample) and upsert_lca_validation
    (SELECT ... FROM lca_validation / INSERT INTO lca_validation)."""

    def __init__(self, sample_species, existing_lca_row, log):
        self.sample_species = sample_species
        self.existing_lca_row = existing_lca_row
        self.log = log
        self._last_query = None

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def execute(self, query, params=None):
        self._last_query = query
        self.log.append((query, params))

    def fetchone(self):
        if "FROM sample" in self._last_query:
            return (self.sample_species,) if self.sample_species is not None else None
        if "FROM lca_validation" in self._last_query:
            return self.existing_lca_row
        return None


class FakeConnection:
    def __init__(self, sample_species, existing_lca_row, log):
        self.cursor_instance = FakeCursor(sample_species, existing_lca_row, log)

    def cursor(self):
        return self.cursor_instance

    def commit(self):
        pass

    def rollback(self):
        pass

    def close(self):
        pass


def make_config(root):
    path = root / "oceanomics.cfg"
    config = configparser.ConfigParser()
    config["postgres"] = {
        "dbname": "test", "user": "u", "password": "p", "host": "h", "port": "5432",
    }
    with path.open("w") as fh:
        config.write(fh)
    return path


class NoNominalSpeciesTests(unittest.TestCase):
    """No nominal_species_id on the sample row: still write the summary file
    (with N/A comparison columns) and still upsert an lca_validation row with
    validated_species_name=None, so PUSH_LCA_BLAST_RESULTS still runs and the
    sample isn't silently absent from lca_validation."""

    def _run(self, force):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prefix = "OG470.ilmn.230607.getorg1770"

            lca_path = root / "lca.12s.tsv"
            with lca_path.open("w", newline="") as fh:
                fh.write("seq_id\tspecies_in_LCA\n")
                fh.write(f"{prefix}.emma102.001\tGenus species\n")

            blast_path = root / "blast.12s.tsv"
            with blast_path.open("w") as fh:
                fh.write(f"{prefix}.emma102\tsubject1\tsome\tblast\tfields\n")

            config_path = make_config(root)

            log = []

            def fake_connect(**_kwargs):
                return FakeConnection(sample_species=None, existing_lca_row=None, log=log)

            class FakePsycopg2:
                connect = staticmethod(fake_connect)

            original_psycopg2 = MODULE.psycopg2
            original_cwd = Path.cwd()
            MODULE.psycopg2 = FakePsycopg2
            os.chdir(root)
            try:
                buffer = io.StringIO()
                with redirect_stdout(buffer):
                    MODULE.compare_lca_and_blast(
                        str(config_path),
                        "OG470",
                        [str(lca_path)],
                        [str(blast_path)],
                        f"lca_results.{prefix}.tsv",
                        assembly_prefix=prefix,
                        force=force,
                    )
            finally:
                MODULE.psycopg2 = original_psycopg2
                os.chdir(original_cwd)

            output = buffer.getvalue()
            with (root / f"lca_results.{prefix}.tsv").open() as fh:
                summary_rows = list(csv.reader(fh, delimiter="\t"))
            return output, summary_rows, log

    def test_writes_na_summary_columns(self):
        output, rows, _log = self._run(force=False)
        self.assertIn("has no nominal_species_id", output)
        self.assertEqual(
            rows[0],
            ["og_id", "LCA_result", "nom_species_id", "Match_YN", "Found_in_blast_YN",
             "validated_rank"],
        )
        self.assertEqual(rows[1][2:], ["N/A", "N/A", "N/A", "N/A"])

    def test_upserts_lca_validation_row_with_null_species(self):
        _output, _rows, log = self._run(force=False)
        insert_calls = [
            (query, params)
            for query, params in log
            if query and "INSERT INTO lca_validation" in query
        ]
        self.assertEqual(len(insert_calls), 1)
        _query, params = insert_calls[0]
        self.assertIsNone(params["validated_species_name"])
        self.assertEqual(params["og_id"], "OG470")
        self.assertEqual(params["tech"], "ilmn")
        self.assertEqual(params["annotation"], "emma102")
        self.assertEqual(params["validator"], "nf-core")


if __name__ == "__main__":
    unittest.main()


class LoadBlastSpeciesSetTests(unittest.TestCase):
    """The BLAST table is parsed, not searched as text.

    The old version returned the whole file lowercased and the caller asked
    whether the nominal name appeared anywhere in it. That failed OPEN, matching
    inside longer names and inside description fields.
    """

    def _write(self, rows):
        import tempfile
        p = Path(tempfile.mkdtemp()) / "blast_combined.tsv"
        p.write_text("".join("\t".join(r) + "\n" for r in rows))
        return str(p)

    def _row(self, sciname, description="desc"):
        # blast_combined is headerless; column 4 (index 3) is the scientific name.
        return ["q", "subj", "1234", sciname, sciname, "Eukaryota", "99.0", description]

    def test_it_collects_species_and_derives_genera(self):
        p = self._write([self._row("Thalassoma lutescens"),
                         self._row("Thalassoma lunare")])
        hits = MODULE.load_blast_species_set(p)
        self.assertEqual(hits.species,
                         frozenset({"thalassoma lutescens", "thalassoma lunare"}))
        self.assertEqual(hits.genera, frozenset({"thalassoma"}))

    def test_self_hits_and_null_names_are_ignored(self):
        p = self._write([self._row("N/A"), self._row("Thalassoma lunare")])
        hits = MODULE.load_blast_species_set(p)
        self.assertEqual(hits.species, frozenset({"thalassoma lunare"}))

    def test_a_name_only_in_a_description_field_does_not_validate(self):
        # The fail-open regression, stated directly. Under the old blob search
        # any string occurring anywhere in the file validated the sample.
        p = self._write([self._row("Thalassoma lunare",
                                   description="voucher of Serrivomer beanii")])
        hits = MODULE.load_blast_species_set(p)
        self.assertNotIn("serrivomer", hits.genera)
        self.assertNotIn("serrivomer beanii", hits.species)

    def test_an_empty_file_yields_empty_sets_without_raising(self):
        hits = MODULE.load_blast_species_set(self._write([]))
        self.assertEqual(hits.species, frozenset())
        self.assertEqual(hits.genera, frozenset())

    def test_a_missing_file_yields_empty_sets_without_raising(self):
        hits = MODULE.load_blast_species_set("/nonexistent/blast_combined.tsv")
        self.assertEqual(hits.species, frozenset())

    def test_an_in_house_reference_hit_is_read_from_its_organism_tag(self):
        # A hit against OceanOmics' own reference database has no NCBI taxid, so
        # columns 3-6 are all 'N/A' and the identification lives only in the
        # subject title as a structured [organism=...] tag. These are real hits,
        # often a sample's only ones -- reading column 4 alone silently discarded
        # nine batch-15 samples that had already reached ENA.
        p = self._write([[
            "q", "NBDL-HK3RDJT1K3246F.v1.mt", "0", "N/A", "N/A", "N/A", "99.5",
            "[organism=Bathypterois oddi] [authority=Sulak, 1977] [mgcode=2] "
            "Bathypterois oddi specimen ANFC H 8091-04 mitochondrion"]])
        hits = MODULE.load_blast_species_set(p)
        self.assertIn("bathypterois oddi", hits.species)
        self.assertIn("bathypterois", hits.genera)

    def test_a_bold_synonym_in_the_subject_title_is_read(self):
        # Taxid 334986 is 'Hydrolagus ogilbyi' to NCBI and 'Chimaera ogilbyi' in
        # BOLD -- the same animal under two combinations. A sample labelled with
        # the synonym is correctly identified, so both names belong in the set.
        p = self._write([[
            "q", "79146_FOAO734-15", "334986", "Hydrolagus ogilbyi",
            "Ogilby's ghostshark", "Eukaryota", "99.7",
            "334986 Chimaera ogilbyi FOAO734-15|Chimaera ogilbyi|COI-5P"]])
        hits = MODULE.load_blast_species_set(p)
        self.assertIn("hydrolagus ogilbyi", hits.species)
        self.assertIn("chimaera ogilbyi", hits.species)

    def test_marker_codes_and_accessions_are_not_mistaken_for_names(self):
        # The binomial shape is what keeps the '|'-delimited parse honest.
        p = self._write([[
            "q", "s", "1", "Squalus mitsukurii", "n", "Eukaryota", "99.0",
            "134996 Squalus mitsukurii FARG335-07|Squalus mitsukurii|COI-5P|EU074610"]])
        hits = MODULE.load_blast_species_set(p)
        self.assertEqual(hits.species, frozenset({"squalus mitsukurii"}))


class MatchAtRankTests(unittest.TestCase):
    """Match at the rank the label asserts, not the loosest one that succeeds."""

    def _hits(self, *names):
        species = frozenset(MODULE.normalise_name(n) for n in names)
        return MODULE.BlastHits(
            species=species,
            genera=frozenset(MODULE.genus_of(s) for s in species))

    def _row(self, genus="dropped", family="dropped"):
        return {"genus": genus, "family": family}

    # -- the positives, one per rank ---------------------------------------

    def test_a_binomial_matches_at_species(self):
        got = MODULE.match_at_rank(("species", "Notacanthus abbotti"),
                                   self._hits("Notacanthus abbotti"), self._row())
        self.assertEqual(got, ("Yes", "species"))

    def test_a_genus_label_matches_when_blast_and_lca_agree(self):
        got = MODULE.match_at_rank(
            ("genus", "Notacanthus"),
            self._hits("Notacanthus abbotti", "Notacanthus bonaparte"),
            self._row(genus="Notacanthus", family="Notacanthidae"))
        self.assertEqual(got, ("Yes", "genus"))

    def test_a_tentative_species_is_released_at_genus_and_recorded_as_downgraded(self):
        got = MODULE.match_at_rank(
            ("species_uncertain", "Squalus"),
            self._hits("Squalus megalops"), self._row(genus="Squalus"))
        self.assertEqual(got, ("Yes", "genus_downgraded"))

    def test_a_family_label_matches_on_the_lca_family(self):
        got = MODULE.match_at_rank(
            ("family", "Ophidiidae"),
            self._hits("Lamprogrammus exutus"),
            self._row(genus="Lamprogrammus", family="Ophidiidae"))
        self.assertEqual(got, ("Yes", "family"))

    # -- the negatives, which are the point --------------------------------

    def test_a_swapped_sample_is_held_at_every_rank(self):
        # THE regression test. A stargazer label whose every BLAST hit and whose
        # LCA both say dogfish shark. The two evidence sources agree with each
        # other and both disagree with the label, so neither can rescue it.
        # If this ever passes, the change is wrong.
        hits = self._hits("Squalus megalops", "Squalus acanthias")
        row = self._row(genus="Squalus", family="Squalidae")
        for nominal in (("species", "Uranoscopus turbisquamatus"),
                        ("genus", "Uranoscopus"),
                        ("species_uncertain", "Uranoscopus"),
                        ("family", "Uranoscopidae")):
            found, rank = MODULE.match_at_rank(nominal, hits, row)
            self.assertEqual(found, "No", nominal)
            self.assertEqual(rank, "unmatched", nominal)

    def test_a_binomial_with_only_genus_support_stays_held(self):
        # This is the design, not a gap. Falling back to genus for a binomial
        # would convert "we could not confirm the species" into "validated".
        got = MODULE.match_at_rank(
            ("species", "Bathypterois parini"),
            self._hits("Bathypterois atricolor"), self._row(genus="Bathypterois"))
        self.assertEqual(got, ("No", "unmatched"))

    def test_a_genus_label_needs_BOTH_blast_and_the_lca(self):
        # BLAST agrees, LCA does not.
        self.assertEqual(
            MODULE.match_at_rank(("genus", "Notacanthus"),
                                 self._hits("Notacanthus abbotti"),
                                 self._row(genus="Squalus")),
            ("No", "unmatched"))
        # LCA agrees, BLAST does not.
        self.assertEqual(
            MODULE.match_at_rank(("genus", "Notacanthus"),
                                 self._hits("Squalus megalops"),
                                 self._row(genus="Notacanthus")),
            ("No", "unmatched"))

    def test_the_lca_absent_sentinels_never_match(self):
        # calculateLCA.py writes 'dropped' where it declined and 'Unknown' where
        # the lineage lookup returned nothing. Neither may be matched as a value.
        for sentinel in ("dropped", "Unknown", "unknown", "", "NA"):
            self.assertEqual(
                MODULE.match_at_rank(("genus", sentinel),
                                     self._hits(f"{sentinel} something"),
                                     self._row(genus=sentinel))[0],
                "No", sentinel)

    def test_an_unparseable_label_never_matches(self):
        self.assertEqual(
            MODULE.match_at_rank((None, None), self._hits("Squalus megalops"),
                                 self._row(genus="Squalus")),
            ("No", "unmatched"))

    def test_empty_evidence_matches_nothing_and_does_not_raise(self):
        empty = MODULE.BlastHits(species=frozenset(), genera=frozenset())
        for nominal in (("species", "Notacanthus abbotti"),
                        ("genus", "Notacanthus"), ("family", "Ophidiidae")):
            self.assertEqual(MODULE.match_at_rank(nominal, empty, {})[0], "No", nominal)

    def test_a_bare_genus_matches_on_the_parsed_genus_set_not_a_substring(self):
        # The fail-open case that reached ENA: 'Serrivomer' used to validate by
        # substring against 'Serrivomer jesperseni'. It still validates -- but at
        # GENUS rank now, recorded as such, and its organism normalises to
        # 'Serrivomer sp.' rather than the bare genus ENA rejects.
        got = MODULE.match_at_rank(
            ("genus", "Serrivomer"), self._hits("Serrivomer jesperseni"),
            self._row(genus="Serrivomer"))
        self.assertEqual(got, ("Yes", "genus"))
        self.assertEqual(
            MODULE.normalise_open_nomenclature("Serrivomer"), "Serrivomer sp.")
