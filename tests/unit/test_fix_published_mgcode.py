"""Unit tests for bin/fix_published_mgcode.py -- the published-tag repair.

The script edits already-published data in place, so the properties that matter
are the conservative ones: it must change the tag and NOTHING else, must never
guess a code it cannot establish, must leave correctly-tagged trees alone, and
must be safe to run twice. Pure stdlib apart from the orf_utils sibling import.
"""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# bin/ scripts import their siblings (orf_utils, ...) the way Nextflow stages
# them: flat on PATH. Mirror that for the file-path load below.
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "fix_published_mgcode", ROOT / "bin" / "fix_published_mgcode.py")
fix = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fix)

GENE_BODY = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTAA"


def gene_fa(tag: int, gene="CO1"):
    return (f">OG1.ilmn.240101.getorg1770|1-39|+|{gene} "
            f"[organism=Acropora millepora] [mgcode={tag}] [topology=linear] "
            f"[gene-coordinates=1-39(+)] Acropora millepora mitochondrially encoded {gene}\n"
            f"{GENE_BODY}\n")


def build_assembly(root: Path, name="OG1.ilmn.240101.getorg1770",
                   true_code=4, tag=2, with_package=True):
    """A minimal published assembly tree: authoritative record + tagged extraction."""
    asm = root / "mitogenomes" / "OG1" / name
    (asm / "annotation").mkdir(parents=True)
    (asm / "genbank" / "processed").mkdir(parents=True)
    (asm / "genbank" / "genes").mkdir(parents=True)
    (asm / "genbank" / "proteins").mkdir(parents=True)
    # authoritative sources: these always carried the right code
    (asm / "genbank" / "processed" / f"{name}.fa").write_text(
        f">{name} [organism=Acropora millepora] [mgcode={true_code}] mitochondrion\n{GENE_BODY}\n")
    (asm / "annotation" / f"{name}.tbl").write_text(
        f"gene\tCO1\n1\t39\tCDS\n\t\t\ttransl_table\t{true_code}\n")
    # extraction outputs: these carried the hardcoded tag
    (asm / "genbank" / "genes" / f"CO1.{name}.fa").write_text(gene_fa(tag))
    (asm / "genbank" / "proteins" / f"CO1.{name}.protein.fa").write_text(gene_fa(tag))
    if with_package:
        pkg = asm / "ena" / "package"
        pkg.mkdir(parents=True)
        (pkg / f"{name}.mitos2110.genes.fa").write_text(gene_fa(tag))
        import hashlib
        h = hashlib.sha256((pkg / f"{name}.mitos2110.genes.fa").read_bytes()).hexdigest()
        (pkg / "checksums.sha256").write_text(f"{h}  {name}.mitos2110.genes.fa\n")
    return asm


class ScanTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_finds_the_mismatched_files(self):
        build_assembly(self.root, true_code=4, tag=2)
        todo, skipped = fix.scan([self.root])
        self.assertEqual(len(todo), 3)          # gene + protein + package copy
        self.assertEqual(skipped, [])
        self.assertTrue(all(new == 4 for _f, _o, new in todo))

    def test_correctly_tagged_tree_is_left_alone(self):
        # A vertebrate: hardcoded 2 was accidentally right, so nothing to do.
        build_assembly(self.root, true_code=2, tag=2)
        todo, _skipped = fix.scan([self.root])
        self.assertEqual(todo, [])

    def test_unestablished_code_is_skipped_not_guessed(self):
        asm = build_assembly(self.root, true_code=4, tag=2)
        # Remove every authoritative source; the tag must NOT be guessed.
        (asm / "genbank" / "processed" / f"{asm.name}.fa").unlink()
        (asm / "annotation" / f"{asm.name}.tbl").unlink()
        todo, skipped = fix.scan([self.root])
        self.assertEqual(todo, [])
        self.assertEqual(len(skipped), 1)

    def test_disagreeing_sources_are_skipped(self):
        asm = build_assembly(self.root, true_code=4, tag=2)
        (asm / "annotation" / f"{asm.name}.tbl").write_text(
            "gene\tCO1\n1\t39\tCDS\n\t\t\ttransl_table\t9\n")   # 4 vs 9
        todo, skipped = fix.scan([self.root])
        self.assertEqual(todo, [])
        self.assertIn("disagree", skipped[0][1])


class RewriteTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_only_the_tag_changes(self):
        asm = build_assembly(self.root, true_code=4, tag=2)
        target = asm / "genbank" / "genes" / f"CO1.{asm.name}.fa"
        before = target.read_text()
        fix.rewrite(target, 4)
        after = target.read_text()
        self.assertIn("[mgcode=4]", after)
        self.assertNotIn("[mgcode=2]", after)
        # everything else byte-identical once the tag is normalised
        norm = lambda s: s.replace("[mgcode=4]", "[mgcode=X]").replace("[mgcode=2]", "[mgcode=X]")
        self.assertEqual(norm(before), norm(after))
        # and the sequence itself is untouched
        self.assertIn(GENE_BODY, after)

    def test_rewrite_is_idempotent(self):
        asm = build_assembly(self.root, true_code=4, tag=2)
        target = asm / "genbank" / "genes" / f"CO1.{asm.name}.fa"
        self.assertTrue(fix.rewrite(target, 4))
        self.assertFalse(fix.rewrite(target, 4))    # second pass is a no-op

    def test_second_scan_after_rewrite_finds_nothing(self):
        build_assembly(self.root, true_code=4, tag=2)
        todo, _ = fix.scan([self.root])
        for f, _o, new in todo:
            fix.rewrite(f, new)
        again, _ = fix.scan([self.root])
        self.assertEqual(again, [])

    def test_checksums_are_refreshed(self):
        import hashlib
        asm = build_assembly(self.root, true_code=4, tag=2)
        pkg = asm / "ena" / "package"
        target = pkg / f"{asm.name}.mitos2110.genes.fa"
        fix.rewrite(target, 4)
        # stale before the refresh
        stale = pkg / "checksums.sha256"
        recorded = stale.read_text().split()[0]
        self.assertNotEqual(recorded, hashlib.sha256(target.read_bytes()).hexdigest())
        fix.refresh_checksums([pkg])
        recorded = stale.read_text().split()[0]
        self.assertEqual(recorded, hashlib.sha256(target.read_bytes()).hexdigest())


class AuthoritativeCodeTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_reads_code_from_the_genome_record(self):
        asm = build_assembly(self.root, true_code=4)
        self.assertEqual(fix.authoritative_code(asm)[0], 4)

    def test_unsupported_code_is_refused(self):
        asm = build_assembly(self.root, true_code=4)
        (asm / "genbank" / "processed" / f"{asm.name}.fa").write_text(
            f">{asm.name} [mgcode=7] mitochondrion\n{GENE_BODY}\n")
        (asm / "annotation" / f"{asm.name}.tbl").write_text(
            "gene\tCO1\n1\t39\tCDS\n\t\t\ttransl_table\t7\n")
        code, note = fix.authoritative_code(asm)
        self.assertIsNone(code)
        self.assertIn("not a supported", note)


if __name__ == "__main__":
    unittest.main()
