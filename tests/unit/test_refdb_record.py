"""Round-trip tests for bin/refdb_record.py.

The per-group <group>_mito_refdb.gb files are no longer tracked -- 84 MB raw /
26 MB compressed across the eight groups, rewritten into git history on every
rebuild, which is why only anthozoa's was tracked and therefore why
sequence-based reference selection was confined to corals. A reference record is
now rebuilt on demand from the four tracked files instead.

These tests are the gate on that substitution: a record rebuilt from
fasta + manifest.tsv + features.tsv must be equivalent, for every field the
pipeline reads, to the same record parsed from the GenBank it was built from.
The fixture is the 3-record slice under
modules/local/select_reference_db/tests/data/, whose NC_083272.1 carries a
2-exon nad5 -- the case that rules out reusing the label database, which stores
each feature's spliced sequence.
"""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# bin/ on sys.path so sibling imports resolve as they do in the task container.
sys.path.insert(0, str(ROOT / "bin"))

FIXTURE = ROOT / "modules" / "local" / "select_reference_db" / "tests" / "data"
SOURCE_GB = FIXTURE / "minidb_source.gb"
REFDB_DIR = FIXTURE / "minidb"

try:
    from Bio import SeqIO  # noqa: F401
    HAVE_BIO = True
except ImportError:
    HAVE_BIO = False


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@unittest.skipUnless(HAVE_BIO, "biopython not installed")
class RefdbRecordRoundTrip(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from Bio import SeqIO
        cls.refdb_record = _load("refdb_record")
        cls.originals = {r.id: r for r in SeqIO.parse(str(SOURCE_GB), "genbank")}
        cls.assertTrue(cls.originals, "fixture GenBank is empty")

    def rebuilt(self, accession):
        return self.refdb_record.materialize(REFDB_DIR, "minidb", accession)

    def test_every_fixture_record_rebuilds(self):
        for acc in self.originals:
            self.assertIsNotNone(self.rebuilt(acc), f"{acc} did not materialise")

    def test_sequence_is_identical(self):
        for acc, orig in self.originals.items():
            self.assertEqual(str(self.rebuilt(acc).seq), str(orig.seq), acc)

    def test_organism_and_lineage_survive(self):
        for acc, orig in self.originals.items():
            got = self.rebuilt(acc)
            self.assertEqual(got.annotations["organism"], orig.annotations["organism"], acc)
            self.assertEqual(got.annotations["taxonomy"], orig.annotations["taxonomy"], acc)

    def test_coral_fix_features_are_identical(self):
        """rrnL, the nad5 exon list and cox1 -- what CORAL_ANNOTATION_FIX transfers."""
        from Bio import SeqIO
        coral_fix_bed = _load("coral_fix_bed")
        for acc, orig in self.originals.items():
            with tempfile.TemporaryDirectory() as td:
                a, b = Path(td) / "orig.gb", Path(td) / "rebuilt.gb"
                SeqIO.write([orig], str(a), "genbank")
                SeqIO.write([self.rebuilt(acc)], str(b), "genbank")
                self.assertEqual(coral_fix_bed.ref_features(a),
                                 coral_fix_bed.ref_features(b), acc)

    def test_divergence_taxonomy_is_identical(self):
        """organism / genus / family / order / lineage, as REFERENCE_DIVERGENCE grades them."""
        from Bio import SeqIO
        rdc = _load("reference_divergence_check")
        for acc, orig in self.originals.items():
            with tempfile.TemporaryDirectory() as td:
                a, b = Path(td) / "orig.gb", Path(td) / "rebuilt.gb"
                SeqIO.write([orig], str(a), "genbank")
                SeqIO.write([self.rebuilt(acc)], str(b), "genbank")
                self.assertEqual(rdc.parse_reference(a), rdc.parse_reference(b), acc)

    def test_multi_exon_feature_keeps_its_exons(self):
        """The reason a coordinates table is needed rather than the label database.

        The label DB stores feat.extract(), i.e. the SPLICED sequence, so an
        intron-split nad5 would come back as one exon and coral_fix_bed's
        exon-wise transfer would silently degrade.
        """
        rebuilt = self.rebuilt("NC_083272.1")
        nad5 = [f for f in rebuilt.features
                if f.type == "CDS"
                and any("ND5" in v.upper() or "NAD5" in v.upper()
                        for v in f.qualifiers.get("gene", []) + f.qualifiers.get("product", []))]
        self.assertTrue(nad5, "no nad5 CDS in the rebuilt record")
        self.assertGreater(len(nad5[0].location.parts), 1,
                           "nad5 came back single-exon; exon structure was lost")


@unittest.skipUnless(HAVE_BIO, "biopython not installed")
class LocationSerialisation(unittest.TestCase):
    def setUp(self):
        self.refdb_record = _load("refdb_record")

    def test_parts_round_trip_including_strand(self):
        for spec in ("0-100:1", "0-100:-1", "10-20:1,30-40:1", "10-20:-1,30-40:-1"):
            loc = self.refdb_record.parse_parts(spec)
            self.assertEqual(self.refdb_record.format_parts(loc), spec)


if __name__ == "__main__":
    unittest.main()
