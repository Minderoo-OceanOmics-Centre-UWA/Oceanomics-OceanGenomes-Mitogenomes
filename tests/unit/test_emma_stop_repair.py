"""EMMA stop repair: extend a PCG onto the stop codon the assembly encodes.

NOVA_260925_LA held 11 Glaucosoma hebraicum because EMMA ended CO1 (9) or ND1 (2)
one codon short of an encoded TAA, calling the last sense codon a poly-A-completed
stop. Four were the same genome as a passing batch-mate at another rotation.

The refusals matter as much as the repairs. EMMA's genuine poly-A stops (a T / TA
tail whose completing bases belong to the next tRNA) are the normal case in a
vertebrate mitogenome and must come back untouched and silent.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))

from orf_utils import revcomp, translate  # noqa: E402

SCRIPT = BIN / "emma_stop_repair.py"
CHROM = "SYN.ilmn.000000.getorg1770"
NAME = f"{CHROM}.emma102"
POLY_A = "putative {} stop codon is completed by the addition of 3' A residues to the mRNA"
BODY = "GCTCTTCCTGTTACT" * 6          # 30 sense codons, no code-2 stop


class Bundle:
    """A tiny EMMA-shaped bundle: genome, .gff, .tbl, cds/, proteins/.

    Features are laid out left to right; each PCG's declared span is written
    exactly as EMMA would (cds/ = span minus its last 3 nt, protein = Met +
    translation from codon 2).
    """

    def __init__(self):
        self.seq = ""
        self.feats = []

    def spacer(self, s):
        self.seq += s

    def pcg(self, gene, coding, after="", strand="+", note=None):
        """coding = the declared span in the gene's own orientation; `after` =
        bases that follow its 3' end in the gene's own orientation."""
        if strand == "+":
            start = len(self.seq) + 1
            self.seq += coding + after
            end = start + len(coding) - 1
        else:
            block = revcomp(coding + after)
            start = len(self.seq) + len(after) + 1
            self.seq += block
            end = start + len(coding) - 1
        self.feats.append(dict(kind="pcg", gene=gene, start=start, end=end,
                               strand=strand, coding=coding, note=note))

    def trna(self, gene, length=20, seq=None):
        start = len(self.seq) + 1
        self.seq += seq if seq else "C" * length
        self.feats.append(dict(kind="trna", gene=gene, start=start,
                               end=len(self.seq), strand="+"))

    def write(self, d):
        d = Path(d)
        (d / "cds").mkdir(parents=True)
        (d / "proteins").mkdir()
        (d / f"{NAME}.fa").write_text(f">{CHROM}\n{self.seq}\n")
        gff = ["##gff-version 3", f"##sequence-region\t{CHROM}\t1\t{len(self.seq)}"]
        tbl = [f">Feature {CHROM}"]
        for i, f in enumerate(self.feats):
            s, e, st, g = f["start"], f["end"], f["strand"], f["gene"]
            a, b = (s, e) if st == "+" else (e, s)
            if f["kind"] == "trna":
                gff.append(f"{CHROM}\tEmma\tgene\t{s}\t{e}\t.\t{st}\t.\tID=g{i};Name={g}")
                gff.append(f"{CHROM}\tEmma\ttRNA\t{s}\t{e}\t1e-20\t{st}\t.\tID=t{i};Parent=g{i};Name={g}")
                tbl += [f"{a}\t{b}\tgene", f"\t\t\tgene\t{g}", f"{a}\t{b}\ttRNA", "\t\t\tproduct\ttRNA-X"]
                continue
            note = f";Note={f['note']}" if f["note"] else ""
            gff.append(f"{CHROM}\tEmma\tgene\t{s}\t{e}\t.\t{st}\t.\tID=g{i};Name={g}")
            gff.append(f"{CHROM}\tEmma\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID=m{i};Parent=g{i};Name={g}")
            gff.append(f"{CHROM}\tEmma\tCDS\t{s}\t{e}\t0.0\t{st}\t0\tID=c{i};Parent=m{i};Name={g};Product=p{note}")
            tbl += [f"{a}\t{b}\tgene", f"\t\t\tgene\t{g}", f"{a}\t{b}\tmRNA", f"\t\t\tgene\t{g}",
                    f"{a}\t{b}\tCDS", "\t\t\tproduct\tp", "\t\t\ttransl_table\t2",
                    f"\t\t\tprotein_id\tgnl|Emma|{i}"]
            if f["note"]:
                tbl.append(f"\t\t\tnote\t{f['note']}")
            cds = f["coding"][:-3]
            (d / "cds" / f"{g}.{NAME}.fa").write_text(f">{NAME}\n{cds}\n")
            (d / "proteins" / f"{g}.{NAME}.fa").write_text(f">{NAME}\nM{translate(cds, 2)[1:]}\n")
        (d / f"{NAME}.gff").write_text("\n".join(gff) + "\n")
        (d / f"{NAME}.tbl").write_text("\n".join(tbl) + "\n")


def run(bundle, dry_run=False):
    tmp = Path(tempfile.mkdtemp())
    bundle.write(tmp / "annotation")
    cmd = [sys.executable, str(SCRIPT), "--annotation-dir", str(tmp / "annotation"),
           "--code", "2", "--status", str(tmp / "s.tsv")]
    if dry_run:
        cmd.append("--dry-run")
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"})
    assert proc.returncode == 0, proc.stderr
    status = [l.split("\t") for l in (tmp / "s.tsv").read_text().splitlines()]
    return tmp / "annotation", status


def cds_row(ann, gene):
    for line in (ann / f"{NAME}.gff").read_text().splitlines():
        f = line.split("\t")
        if len(f) > 8 and f[2] == "CDS" and f"Name={gene};" in f[8] + ";":
            return int(f[3]), int(f[4]), f[8]
    raise AssertionError(gene)


def fasta(path):
    return "".join(l.strip() for l in path.read_text().splitlines() if not l.startswith(">"))


class TestRepairs(unittest.TestCase):

    def test_sense_codon_called_stop_is_extended_onto_encoded_taa(self):
        # The CO1 signature: "...GCA AAA|TAA", span % 3 == 0, note claims AAA.
        b = Bundle()
        b.spacer("GG")
        b.pcg("MT-CO1", "ATG" + BODY + "GCAAAA", after="TAA", note=POLY_A.format("AAA"))
        b.spacer("GGGGG")
        b.trna("MT-TS2")
        ann, status = run(b)
        self.assertEqual(status, [["REPAIRED", "CO1",
                                   "stop_extended(+3,TAA) 3-101 -> 3-104 +"]])
        s, e, attrs = cds_row(ann, "MT-CO1")
        self.assertEqual((s, e), (3, 104))
        self.assertNotIn("Note=", attrs)
        tbl = (ann / f"{NAME}.tbl").read_text()
        self.assertEqual(tbl.count("3\t104\t"), 3)        # gene, mRNA, CDS
        self.assertNotIn("completed by", tbl)
        self.assertEqual(fasta(ann / "cds" / f"MT-CO1.{NAME}.fa"), "ATG" + BODY + "GCAAAA")
        self.assertTrue(fasta(ann / "proteins" / f"MT-CO1.{NAME}.fa").endswith("AK"))

    def test_partial_ta_completed_by_encoded_a_is_extended_one_base(self):
        # The ND1 974 bp signature (OG3174, OG3237): "...CTC TA|A", gap before TI.
        b = Bundle()
        b.pcg("MT-ND1", "ATG" + BODY + "CTCTA", after="A", note=POLY_A.format("TAA"))
        b.spacer("GGGG")
        b.trna("MT-TI")
        ann, status = run(b)
        self.assertEqual(status[0][:2], ["REPAIRED", "ND1"])
        self.assertIn("stop_extended(+1,TAA)", status[0][2])
        self.assertEqual(fasta(ann / "cds" / f"MT-ND1.{NAME}.fa"), "ATG" + BODY + "CTC")

    def test_minus_strand_gene_extends_its_low_coordinate(self):
        b = Bundle()
        b.trna("MT-TE")
        b.spacer("GGGG")
        b.pcg("MT-ND6", "ATG" + BODY + "CCCAAA", after="TAG", strand="-",
              note=POLY_A.format("AAA"))
        ann, _ = run(b)
        s, e, _attrs = cds_row(ann, "MT-ND6")
        old_s = 20 + 4 + 3 + 1
        self.assertEqual(s, old_s - 3)
        tbl = (ann / f"{NAME}.tbl").read_text()
        self.assertIn(f"{e}\t{s}\tCDS", tbl)               # tbl writes minus as high->low


class TestRefusals(unittest.TestCase):

    def test_genuine_poly_a_stop_into_next_trna_is_untouched_and_silent(self):
        # COX2 ending in T with tRNA-Lys starting "AA": the TAA is not encoded
        # outside another feature, so EMMA's poly-A stop is right.
        b = Bundle()
        b.pcg("MT-CO2", "ATG" + BODY + "T", note=POLY_A.format("TAA"))
        b.trna("MT-TK", seq="AA" + "C" * 18)
        before = Path(tempfile.mkdtemp())
        b.write(before / "a")
        ann, status = run(b)
        self.assertEqual(status, [["NONE", "-", "no truncated stops"]])
        self.assertEqual((ann / f"{NAME}.gff").read_text(),
                         (before / "a" / f"{NAME}.gff").read_text())

    def test_encoded_stop_inside_neighbouring_feature_is_refused(self):
        b = Bundle()
        b.pcg("MT-ND3", "ATG" + BODY + "GCAAAA", note=POLY_A.format("AAA"))
        b.trna("MT-TR", seq="TAA" + "C" * 17)
        ann, status = run(b)
        self.assertEqual(status, [["SKIP", "ND3", "overlaps_MT-TR"]])
        self.assertEqual(cds_row(ann, "MT-ND3")[:2], (1, 99))

    def test_aga_agg_are_not_accepted_as_a_new_stop(self):
        b = Bundle()
        b.pcg("MT-CO1", "ATG" + BODY + "GCAAAA", after="AGG", note=POLY_A.format("AAA"))
        b.spacer("GGGGG")
        ann, status = run(b)
        self.assertEqual(status, [["SKIP", "CO1", "no_encoded_stop(AGG)"]])

    def test_complete_stop_is_left_alone(self):
        b = Bundle()
        b.pcg("MT-ND5", "ATG" + BODY + "TAA")
        b.spacer("GGGG")
        _ann, status = run(b)
        self.assertEqual(status, [["NONE", "-", "no truncated stops"]])

    def test_extension_past_linearisation_point_is_refused(self):
        b = Bundle()
        b.pcg("MT-CYTB", "ATG" + BODY + "GCAAAA", note=POLY_A.format("AAA"))
        _ann, status = run(b)
        self.assertEqual(status, [["SKIP", "CYTB", "extension_crosses_origin"]])

    def test_dry_run_reports_without_writing(self):
        b = Bundle()
        b.pcg("MT-CO1", "ATG" + BODY + "GCAAAA", after="TAA", note=POLY_A.format("AAA"))
        b.spacer("GGG")
        ann, status = run(b, dry_run=True)
        self.assertEqual(status[0][0], "REPAIRED")
        self.assertEqual(cds_row(ann, "MT-CO1")[:2], (1, 99))

    def test_incomplete_bundle_is_skipped_not_crashed(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "annotation").mkdir()
        proc = subprocess.run([sys.executable, str(SCRIPT), "--annotation-dir",
                               str(tmp / "annotation"), "--code", "2"],
                              capture_output=True, text=True,
                              env={"PYTHONPATH": str(BIN), "PATH": "/usr/bin:/bin"})
        self.assertEqual(proc.returncode, 0)
        self.assertTrue(proc.stdout.startswith("SKIP"))


if __name__ == "__main__":
    unittest.main()
