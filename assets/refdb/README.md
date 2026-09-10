# Invertebrate mitogenome reference databases

One curated database per taxon group, built by
[`bin/build_invert_reference_db.py`](../../bin/build_invert_reference_db.py) from NCBI
RefSeq complete mitogenomes.

Resolving what a sample is seeded from, or referenced against, is **two narrowing
stages**. Both consumers below share stage 2:

1. **By taxonomy** — `InvertTaxonGroups.seedDbGroup()` in
   [`lib/`](../../lib/InvertTaxonGroups.groovy) maps the sample's class to one of the
   groups here. A class in no group is **not** reseeded and keeps its first-pass
   assembly, rather than being seeded from a wrong phylum.
2. **By sequence** — `SELECT_REFERENCE_DB` BLASTs the sample's own assembly against
   `<group>_mito_refdb.fasta` and ranks the records by total bitscore. Label-free, so a
   wrong or coarse species label cannot decide the answer.

The two consumers:

- **GetOrganelle reseed** (all groups) — `GETORGANELLE_RESEED` seeds from the
  **top `params.reseed_seed_top_n`** records (`-s`) and labels contigs with just those
  records' genes (`--genes`). It used to be handed the whole group, which recruits reads
  from across the phylum, halves effective coverage and shatters the graph: that took
  INV04_BOLOCERA's reseed from 2 scaffolds to 12. Top-n rather than the single best,
  because from a small fragmented first pass the pick is reliable at order/subclass level
  but not at species. A sample nothing aligns to is not reseeded — that is the honest
  signal that the group does not represent its lineage, and falling back to the whole
  group is the very failure this avoids.
- **Annotation reference** (all groups) — the single best record becomes the reference
  for `CORAL_ANNOTATION_FIX`, `GETORGANELLE_CHECK` and `REFERENCE_RELEVANCE`.
  `anthozoa_reference.gb` is the single curated fallback for the rare coral no DB record
  aligns to. Note that a reference for every group does not put every group through the
  coral fixer: that stays gated on `CORAL_FIX_ELIGIBLE_CLASSES` (Cnidaria), because the
  nad5-717 intron repair it performs is Hexacorallia-specific.

## Contents

| group | Entrez organism expression | min CDS | needs both rRNAs | needs nad5 | RefSeq only | records | families |
|---|---|---|---|---|---|---|---|
| anthozoa | `txid6101` | 13 | yes | yes | yes | 221 | 87 |
| porifera | `txid6040` | 13 | yes | yes | yes | 58 | 33 |
| mollusca | `txid6447` | 12 | yes | no | yes | 850 | 208 |
| arthropoda | `txid6657` NOT Hexapoda/Arachnida/Myriapoda | 13 | yes | no | yes | 647 | 169 |
| echinodermata | `txid7586` | 13 | yes | no | yes | 135 | 56 |
| ctenophora | `txid10197` | 10 | no | no | **no** | **16** | **7** |
| tunicata | `txid7712` | 12 | yes | no | yes | 19 | 7 |
| annelida | `txid6340` | 12 | yes | no | yes | 167 | 50 |

Built 2026-09-02, except **anthozoa**, which is the original 2026-08 build kept unchanged:
rebuilding it today yields 278 records (a strict superset of the 221 — no record is lost,
57 are added by new RefSeq entries and the widened rRNA-synonym matcher), but those extra
records also change which reference `SELECT_REFERENCE_DB` picks for every coral, so
that refresh belongs in its own change with its own coral annotation check.

The `.features.tsv` and `lineage` column were added later, with
`--refresh-derived` (below) rather than a rebuild, precisely so that schema change could
not smuggle in the content change above. `.fasta` and `.label.fasta` came back
byte-identical for all eight groups, and all 2,101 records round-trip through
`refdb_record.py` to identical `coral_fix_bed.ref_features()` and
`reference_divergence_check.parse_reference()` output.

The completeness bar is per group and not negotiable upward for its own sake: coral
references must carry the features `CORAL_ANNOTATION_FIX` transfers (both rRNAs, a nad5
CDS) plus the 13-PCG cnidarian set, while ctenophore mitogenomes are genuinely reduced
(no atp6, no tRNAs, ~10 PCGs, rRNAs often unannotated), so the coral bar would reject
every valid ctenophore record.

## The RefSeq restriction, and why ctenophora is exempt

Every group but ctenophora searches `AND refseq[filter]`. RefSeq is a curated
subset, not a completeness bar — the records it omits are ordinary INSDC
submissions that pass exactly the same `record_is_complete()` check. For a
well-populated group the restriction is a useful de-duplicator and stays on.

For ctenophora it was the binding constraint, and it hid most of the phylum:
4 records out of 35 matching, with **no Platyctenida at all**. That is what made
`INV08_TJALFIELLA` look like an unfixable reference gap in the 20-sample invert
panel — it was seeded from three families in two other orders while **two
*Tjalfiella* mitogenomes**, its own genus, sat in GenBank behind the filter
(PP327218, 11,397 bp, 11 CDS; PP331237, 11,020 bp, 11 CDS — both clear the
group's own min-CDS bar of 10). Lifting it takes the group 4 → 16 records and
3 → 7 families, adding Tjalfiellidae, Lyroctenidae, Benthoplanidae and
Euplokamidae.

It is deliberately per group, not global. Lifting it everywhere widens the other
seven 2.5–4.7× (mollusca 855 → 3169, anthozoa 295 → 1280), which re-picks the
`SELECT_REFERENCE_DB` reference for samples that are already submitted — the same
risk that keeps the anthozoa build frozen above. Flip one with `refseq_only` in
`GROUPS`, and only with its own before/after reference-selection diff.

Two build-time dedup stages exist only because of this, and are no-ops while the
filter is on:

- **INSDC twins.** Without the filter a RefSeq record arrives alongside the
  submission it was derived from (`NC_038065` + `MG655622`). `drop_insdc_twins()`
  reads the `reference sequence is identical to X` comment and keeps the RefSeq copy.
- **Per-organism cap** (`--max-per-organism`, default 2). The widened ctenophore
  search returns nine *Vallicula multiformis* isolates. All nine sit in one family,
  so neither a top-n seed panel nor `select_fallback_seed`'s family balancing can
  dilute them. Ranked RefSeq first, then longest, then most CDS, so it is reproducible.

Arthropoda subtracts Hexapoda (`txid6960`), Arachnida (`txid6854`) and Myriapoda
(`txid61985`). Those are >95% of arthropod RefSeq mitogenomes, none of them is anything
OceanOmics sequences, and leaving them in both truncated the search at `--retmax` and
buried the marine crustaceans and pycnogonids the seed exists for.

## What is tracked

Four files per group: `.fasta`, `.label.fasta`, `.manifest.tsv` and `.features.tsv`.

The `.gb` files are **not** tracked for any group. They are 84 MB raw / 26 MB compressed
across the eight, they are plain git blobs (no LFS), and every rebuild rewrites all of
them into history. Everything the pipeline reads out of a reference GenBank —  sequence
and length, organism and taxonomy, and gene features with their exon structure — is in
the four tracked files, and [`bin/refdb_record.py`](../../bin/refdb_record.py) rebuilds
an equivalent record on demand. Downstream consumers are unchanged: they still take a
`--ref-gb` / `--reference-gb` path, which also keeps them working for vertebrates, whose
reference really is a GenBank downloaded by `findMitoReference`.

`.features.tsv` exists rather than reusing `.label.fasta` because the label database
stores each feature's **spliced** sequence (`feat.extract`), while `coral_fix_bed.py`
needs one entry per exon of the group-I-intron-split nad5. Coordinates keep the exon
structure; the sequence comes from the `.fasta`.

The `.nucl.*` BLAST databases are gitignored too — see [`.gitignore`](.gitignore).
Regenerate any untracked file by rerunning the builder.

## Rebuilding

Needs biopython (and `makeblastdb` on PATH for the BLAST db). On Setonix there is a venv at
`/software/projects/pawsey1348/tpeirce/refdb_venv`:

```bash
PY=/software/projects/pawsey1348/tpeirce/refdb_venv/bin/python
$PY bin/build_invert_reference_db.py --group mollusca --email you@uwa.edu.au   # one group
$PY bin/build_invert_reference_db.py --all --email you@uwa.edu.au              # all of them
```

To regenerate the tracked files from a group's existing `.gb` **without** contacting NCBI
— a schema change, say — use `--refresh-derived`. A fresh download would also change which
records the database holds, so this is the way to refresh the artifacts alone:

```bash
$PY bin/build_invert_reference_db.py --all --refresh-derived
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
3. Check the four tracked files are all present (`.gb` should not be committed).
4. Update the table above.

In-house (non-RefSeq) assemblies are not included in these shareable builds. Merge them
for a private build with `--extra-gb <files-or-dirs>`; they are held to the same bar.
