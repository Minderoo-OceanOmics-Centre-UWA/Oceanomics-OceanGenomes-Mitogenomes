# Invertebrate mitogenome reference databases

One curated database per taxon group, built by
[`bin/build_invert_reference_db.py`](../../bin/build_invert_reference_db.py) from NCBI
RefSeq complete mitogenomes. Two consumers:

- **GetOrganelle reseed** (all groups) — when an invertebrate's first-pass assembly
  fails, `GETORGANELLE_RESEED` seeds from `<group>_mito_refdb.fasta` (`-s`) and labels
  contigs with `<group>_mito_refdb.label.fasta` (`--genes`). The group is resolved from
  the sample's class by `InvertTaxonGroups.seedDbGroup()` in [`lib/`](../../lib/InvertTaxonGroups.groovy);
  a class in no group is **not** reseeded and keeps its first-pass assembly, rather than
  being seeded from a wrong phylum.
- **Annotation reference** (anthozoa only, so far) — `SELECT_CORAL_REFERENCE` BLASTs the
  assembly against `anthozoa_mito_refdb.gb` and picks a reference by sequence similarity,
  so a wrong or coarse species label can no longer hand the coral fixer a wrong-family
  reference. `anthozoa_reference.gb` is the single curated fallback for the rare sample
  no DB record aligns to.

## Contents

| group | Entrez organism expression | min CDS | needs both rRNAs | needs nad5 | records | families |
|---|---|---|---|---|---|---|
| anthozoa | `txid6101` | 13 | yes | yes | 221 | 87 |
| porifera | `txid6040` | 13 | yes | yes | 58 | 33 |
| mollusca | `txid6447` | 12 | yes | no | 850 | 208 |
| arthropoda | `txid6657` NOT Hexapoda/Arachnida/Myriapoda | 13 | yes | no | 647 | 169 |
| echinodermata | `txid7586` | 13 | yes | no | 135 | 56 |
| ctenophora | `txid10197` | 10 | no | no | 4 | 3 |
| tunicata | `txid7712` | 12 | yes | no | 19 | 7 |
| annelida | `txid6340` | 12 | yes | no | 167 | 50 |

Built 2026-09-02, except **anthozoa**, which is the original 2026-08 build kept unchanged:
rebuilding it today yields 278 records (a strict superset of the 221 — no record is lost,
57 are added by new RefSeq entries and the widened rRNA-synonym matcher), but those extra
records also change which reference `SELECT_CORAL_REFERENCE` picks for every coral, so
that refresh belongs in its own change with its own coral annotation check.

The completeness bar is per group and not negotiable upward for its own sake: coral
references must carry the features `CORAL_ANNOTATION_FIX` transfers (both rRNAs, a nad5
CDS) plus the 13-PCG cnidarian set, while ctenophore mitogenomes are genuinely reduced
(no atp6, no tRNAs, ~10 PCGs, rRNAs often unannotated), so the coral bar would reject
every valid ctenophore record.

Arthropoda subtracts Hexapoda (`txid6960`), Arachnida (`txid6854`) and Myriapoda
(`txid61985`). Those are >95% of arthropod RefSeq mitogenomes, none of them is anything
OceanOmics sequences, and leaving them in both truncated the search at `--retmax` and
buried the marine crustaceans and pycnogonids the seed exists for.

## What is tracked

`.fasta`, `.label.fasta` and `.manifest.tsv` for every group; `.gb` for anthozoa only.
The other `.gb` files (85 MB across the groups, read by nothing today) and the
`.nucl.*` BLAST databases are gitignored — see [`.gitignore`](.gitignore). Regenerate any
of them by rerunning the builder.

## Rebuilding

Needs biopython (and `makeblastdb` on PATH for the BLAST db). On Setonix there is a venv at
`/software/projects/pawsey1348/tpeirce/refdb_venv`:

```bash
PY=/software/projects/pawsey1348/tpeirce/refdb_venv/bin/python
$PY bin/build_invert_reference_db.py --group mollusca --email you@uwa.edu.au   # one group
$PY bin/build_invert_reference_db.py --all --email you@uwa.edu.au              # all of them
```

Each group writes to `assets/refdb/<group>/` by default. The builder aborts if the Entrez
search returns more records than `--retmax`, so a group can never be silently built from
an arbitrary slice of itself.

## Adding a group

1. Add it to `GROUPS` in `bin/build_invert_reference_db.py` with its organism expression
   and completeness bar, and build it.
2. Add its classes to `InvertTaxonGroups` and return the group from `seedDbGroup()` —
   the group key and the directory name must match, or the reseed looks for a database
   that does not exist.
3. Update the table above.

In-house (non-RefSeq) assemblies are not included in these shareable builds. Merge them
for a private build with `--extra-gb <files-or-dirs>`; they are held to the same bar.
