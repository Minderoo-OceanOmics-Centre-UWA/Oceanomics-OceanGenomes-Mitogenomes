"""Unit tests for bin/mito_depth.py.

Driven entirely through --sam with synthetic SAM records, so no minimap2 is
needed and the tests run anywhere. Covers the three things most likely to be
subtly wrong:

  * the circular fold (must reconstruct a uniform profile, and must NOT double
    count an origin-spanning read),
  * gap-compressed identity (a real control-region indel must survive a filter
    that a raw NM/aligned-length ratio would fail it on),
  * the guards that decide when doubling is legitimate at all.
"""
import csv
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "bin" / "mito_depth.py"


def genome(n=100, seed=3):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def sam_record(name, pos, cigar, seq_len, ref="chrM", flag=0, nm=None, de=None):
    """One minimal but valid SAM line."""
    fields = [
        name, str(flag), ref, str(pos), "60", cigar, "*", "0", "0",
        "A" * seq_len, "I" * seq_len,
    ]
    if nm is not None:
        fields.append("NM:i:{}".format(nm))
    if de is not None:
        fields.append("de:f:{}".format(de))
    return "\t".join(fields)


class MitoDepthTests(unittest.TestCase):
    def run_depth(self, seq, sam_lines, circular="true", fasta_name="OG1.hifi.d.asm.fasta",
                  min_identity=0.95, min_aligned_frac=0.80, records=None):
        """Run mito_depth.py over a synthetic SAM and return the parsed row."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = tmp / fasta_name
            with fasta.open("w") as handle:
                if records is None:
                    handle.write(">chrM\n")
                    for i in range(0, len(seq), 60):
                        handle.write(seq[i:i + 60] + "\n")
                else:
                    for name, sub in records:
                        handle.write(">{}\n".format(name))
                        for i in range(0, len(sub), 60):
                            handle.write(sub[i:i + 60] + "\n")

            sam = tmp / "in.sam"
            sam.write_text("@SQ\tSN:chrM\tLN:{}\n".format(len(seq)) + "\n".join(sam_lines) + "\n")

            out = tmp / "out.mito_depth.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPT),
                 "--fasta", str(fasta), "--sample", "OG1.hifi.d.asm",
                 "--out", str(out), "--sam", str(sam),
                 "--circular", circular, "--workdir", str(tmp),
                 "--min-identity", str(min_identity),
                 "--min-aligned-frac", str(min_aligned_frac)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with out.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(rows, "no data row written; stderr={}".format(result.stderr))
            return rows[0]

    # --- circular fold -----------------------------------------------------

    def test_fold_gives_uniform_depth(self):
        """Two reads tiling a 100 bp circle end to end -> flat depth 1x."""
        seq = genome(100)
        # On the doubled (200 bp) reference: 1-50 and 151-200 (i.e. 51-100 folded).
        sam = [
            sam_record("a", 1, "50M", 50, nm=0),
            sam_record("b", 151, "50M", 50, nm=0),
        ]
        row = self.run_depth(seq, sam)
        self.assertEqual(row["circular_doubled"], "true")
        self.assertEqual(row["mean_depth"], "1")
        self.assertEqual(row["depth_cv"], "0")
        self.assertEqual(row["breadth_1x"], "1")
        self.assertEqual(row["target_length_bp"], "100")

    def test_junction_spanning_read_not_double_counted(self):
        """A read across the origin must add exactly 1x, not 2x, anywhere."""
        seq = genome(100)
        # 91..110 on the doubled reference = 91..100 plus 1..10 of the monomer.
        row = self.run_depth(seq, [sam_record("j", 91, "20M", 20, nm=0)])
        self.assertEqual(row["circular_doubled"], "true")
        # 20 covered bases over a 100 bp monomer, none of them stacked.
        self.assertEqual(row["mean_depth"], "0.2")
        self.assertEqual(row["breadth_1x"], "0.2")
        self.assertEqual(row["min_depth"], "0")
        # Nothing anywhere may reach 2x.
        self.assertEqual(row["breadth_10x"], "0")

    # --- identity filtering ------------------------------------------------

    def test_gap_compressed_identity_keeps_indel_read(self):
        """75M30D45M with NM=30: raw identity 0.75, gap-compressed 0.99.

        A raw NM/aligned-length filter at 0.95 would discard this real
        control-region indel read. It must be kept.
        """
        seq = genome(400)
        row = self.run_depth(
            seq, [sam_record("d", 1, "75M30D45M", 120, nm=30)],
            circular="false", min_aligned_frac=0.0,
        )
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["reads_fail_identity"], "0")

    def test_low_identity_read_excluded(self):
        seq = genome(400)
        sam = [sam_record("numt", 1, "150M", 150, nm=15)]  # identity 0.90
        row = self.run_depth(seq, sam, circular="false", min_aligned_frac=0.0)
        self.assertEqual(row["mito_mapped_reads"], "0")
        self.assertEqual(row["reads_fail_identity"], "1")

    def test_low_identity_read_kept_at_looser_threshold(self):
        seq = genome(400)
        sam = [sam_record("numt", 1, "150M", 150, nm=15)]
        row = self.run_depth(seq, sam, circular="false", min_identity=0.85,
                             min_aligned_frac=0.0)
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["reads_fail_identity"], "0")

    def test_de_tag_preferred_over_nm(self):
        """de:f: wins when both are present: NM alone would fail this read."""
        seq = genome(400)
        sam = [sam_record("x", 1, "150M", 150, nm=45, de=0.001)]
        row = self.run_depth(seq, sam, circular="false", min_aligned_frac=0.0)
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["reads_fail_identity"], "0")

    # --- clipping ----------------------------------------------------------

    def test_soft_clipped_read_excluded_in_middle(self):
        seq = genome(400)
        # 60 of 150 bases aligned, well inside the reference -> chimeric.
        sam = [sam_record("clip", 100, "60M90S", 150, nm=0)]
        row = self.run_depth(seq, sam, circular="false")
        self.assertEqual(row["mito_mapped_reads"], "0")
        self.assertEqual(row["reads_fail_clip"], "1")

    def test_clip_waived_at_reference_end(self):
        """Same clip fraction, but the reference ran out -> keep it."""
        seq = genome(400)
        sam = [sam_record("edge", 341, "60M90S", 150, nm=0)]
        row = self.run_depth(seq, sam, circular="false")
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["reads_fail_clip"], "0")

    # --- doubling guards ---------------------------------------------------

    def test_no_doubling_when_not_circular(self):
        seq = genome(100)
        row = self.run_depth(seq, [sam_record("a", 1, "50M", 50, nm=0)], circular="false")
        self.assertEqual(row["circular_doubled"], "false")

    def test_no_doubling_for_concat_assembly(self):
        """A _concat molecule is a blind join of scaffolds, never a circle."""
        seq = genome(100)
        row = self.run_depth(seq, [sam_record("a", 1, "50M", 50, nm=0)],
                             circular="true", fasta_name="OG1.hifi.d.asm_concat.fasta")
        self.assertEqual(row["circular_doubled"], "false")

    def test_no_doubling_for_multi_record_fasta(self):
        seq = genome(100)
        records = [("c1", seq[:50]), ("c2", seq[50:])]
        row = self.run_depth(seq, [sam_record("a", 1, "20M", 20, nm=0)],
                             circular="true", records=records)
        self.assertEqual(row["circular_doubled"], "false")
        self.assertEqual(row["n_contigs"], "2")

    # --- bookkeeping -------------------------------------------------------

    def test_unmapped_counted_in_total_but_not_depth(self):
        seq = genome(400)
        sam = [
            sam_record("m", 1, "100M", 100, nm=0),
            sam_record("u", 0, "*", 100, flag=4),
        ]
        row = self.run_depth(seq, sam, circular="false", min_aligned_frac=0.0)
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["total_reads"], "2")
        self.assertEqual(row["mito_read_fraction"], "0.5")

    def test_secondary_and_supplementary_ignored(self):
        seq = genome(400)
        sam = [
            sam_record("m", 1, "100M", 100, nm=0),
            sam_record("m", 200, "100M", 100, flag=256, nm=0),
            sam_record("m", 300, "50M", 50, flag=2048, nm=0),
        ]
        row = self.run_depth(seq, sam, circular="false", min_aligned_frac=0.0)
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["total_reads"], "1")
        self.assertEqual(row["supplementary_dropped"], "1")

    def test_deletion_does_not_contribute_depth(self):
        """samtools depth semantics: a deleted base is not covered."""
        seq = genome(400)
        row = self.run_depth(seq, [sam_record("d", 1, "50M50D50M", 100, nm=50)],
                             circular="false", min_aligned_frac=0.0)
        # 100 covered bases of 400, the 50 deleted ones excluded.
        self.assertEqual(row["breadth_1x"], "0.25")

    def test_depth_method_and_alias_columns(self):
        seq = genome(100)
        row = self.run_depth(seq, [sam_record("a", 1, "50M", 50, nm=0)], circular="false")
        self.assertEqual(row["depth_method"], "remap_full_v1")
        # The summary reads these aliases; they must track the depth columns.
        self.assertEqual(row["mean_coverage"], row["mean_depth"])
        self.assertEqual(row["coverage_cv"], row["depth_cv"])

    def test_column_names_do_not_trip_numt_detector(self):
        """has_numt_signal() greps every .tsv for these substrings."""
        forbidden = ["rejected", "divergent", "low_coverage", "low confidence",
                     "nuclear mitochondrial"]
        header = "\t".join(_load_columns()).lower()
        for pattern in forbidden:
            self.assertNotIn(pattern, header)

    # --- regressions found in review ---------------------------------------

    def test_insertion_counts_toward_aligned_fraction(self):
        """50M100I50M is fully aligned; only the reference lacks the insert.

        Counting I as unaligned made it look 50% clipped, so a read carrying a
        large insertion failed the fraction floor -- the very case gap-compressed
        identity exists to keep.
        """
        seq = genome(400)
        row = self.run_depth(
            seq, [sam_record("ins", 100, "50M100I50M", 200, nm=100, de=0.005)],
            circular="false", min_aligned_frac=0.80,
        )
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["reads_fail_clip"], "0")

    def test_edge_waiver_is_directional(self):
        """1M149S at position 1 must NOT be waived: the clip faces inward."""
        seq = genome(400)
        row = self.run_depth(seq, [sam_record("junk", 1, "1M149S", 150, nm=0)],
                             circular="false", min_aligned_frac=0.80)
        self.assertEqual(row["mito_mapped_reads"], "0")
        self.assertEqual(row["reads_fail_clip"], "1")

    def test_edge_waiver_accepts_clip_running_off_the_start(self):
        """90S60M at position 1: the leading clip really did run off the start."""
        seq = genome(400)
        row = self.run_depth(seq, [sam_record("edge", 1, "90S60M", 150, nm=0)],
                             circular="false", min_aligned_frac=0.80)
        self.assertEqual(row["mito_mapped_reads"], "1")
        self.assertEqual(row["reads_fail_clip"], "0")

    def test_mapper_failure_does_not_report_zero_depth(self):
        """A broken mapper must fail open, never look like a real 0x measurement.

        Writing mean_depth=0 with depth_method=remap_full_v1 would be
        indistinguishable from a genuine result and would reach the database.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = tmp / "asm.fasta"
            fasta.write_text(">chrM\n" + genome(200) + "\n")
            reads = tmp / "r.fastq"
            reads.write_text("@a\nACGT\n+\nIIII\n")
            out = tmp / "out.mito_depth.tsv"
            env = dict(os.environ)
            env["PATH"] = str(tmp)  # no minimap2, no awk on PATH
            result = subprocess.run(
                [sys.executable, str(SCRIPT),
                 "--fasta", str(fasta), "--reads", str(reads),
                 "--sample", "OG1", "--out", str(out),
                 "--preset", "sr", "--circular", "false", "--workdir", str(tmp)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True, env=env,
            )
            self.assertEqual(result.returncode, 0)
            lines = out.read_text().strip().split("\n")
            self.assertEqual(len(lines), 1, "expected header only, got: {}".format(lines))

    def test_failure_writes_header_only_and_exits_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "out.mito_depth.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPT),
                 "--fasta", str(tmp / "missing.fasta"), "--sample", "OG1",
                 "--out", str(out), "--sam", str(tmp / "also_missing.sam"),
                 "--workdir", str(tmp)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )
            self.assertEqual(result.returncode, 0)
            self.assertTrue(out.exists())
            lines = out.read_text().strip().split("\n")
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith("sample\t"))


class SubsampleTests(unittest.TestCase):
    """The hash draw must be uniform and keep mates together."""

    def test_hash_draw_keeps_mates_together_and_hits_rate(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1 = tmp / "r1.fastq"
            r2 = tmp / "r2.fastq"
            n = 4000
            for path, mate in ((r1, "1"), (r2, "2")):
                with path.open("w") as handle:
                    for i in range(n):
                        handle.write("@read{}/{}\nACGT\n+\nIIII\n".format(i, mate))

            out1 = tmp / "s1.fastq"
            out2 = tmp / "s2.fastq"
            kept1 = module.subsample_fastq([str(r1)], str(out1), 0.25)
            kept2 = module.subsample_fastq([str(r2)], str(out2), 0.25)

            # Same reads chosen from both mates.
            self.assertEqual(kept1, kept2)
            names1 = [l[1:].strip().rsplit("/", 1)[0]
                      for l in out1.read_text().split("\n") if l.startswith("@read")]
            names2 = [l[1:].strip().rsplit("/", 1)[0]
                      for l in out2.read_text().split("\n") if l.startswith("@read")]
            self.assertEqual(names1, names2)

            # Uniform enough: 25% of 4000 within a generous tolerance.
            self.assertGreater(kept1, n * 0.20)
            self.assertLess(kept1, n * 0.30)


def _load_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("mito_depth", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_columns():
    return _load_module().COLUMNS


if __name__ == "__main__":
    unittest.main()
