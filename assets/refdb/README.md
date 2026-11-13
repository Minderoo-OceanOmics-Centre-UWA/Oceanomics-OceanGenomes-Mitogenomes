# Invertebrate mitogenome reference databases

One curated database per taxon group, built by
[`bin/build_invert_reference_db.py`](../../bin/build_invert_reference_db.py) from NCBI
complete mitogenomes — all of INSDC, not only the RefSeq subset.

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

| group | Entrez organism expression | min CDS | needs both rRNAs | needs nad5 | records | families | genera |
|---|---|---|---|---|---|---|---|
| anthozoa | `txid6101` | 13 | yes | yes | 638 | 124 | 308 |
| porifera | `txid6040` | 13 | yes | yes | 87 | 41 | 60 |
| mollusca | `txid6447` | 12 | yes | no | 1391 | 258 | 767 |
| arthropoda | `txid6657` NOT Hexapoda/Arachnida/Myriapoda | 13 | yes | no | 1078 | 231 | 561 |
| echinodermata | `txid7586` | 13 | yes | no | 263 | 80 | 171 |
| ctenophora | `txid10197` | 10 | no | no | 16 | 7 | 9 |
| tunicata | `txid7712` | 12 | yes | no | 33 | 11 | 21 |
| annelida | `txid6340` | 12 | yes | no | 420 | 63 | 208 |

Every group searches all of INSDC and keeps at most one record per organism
(ctenophora keeps two — see below). Rebuilt 2026-09-10 from the RefSeq-only builds
below, which is where the `refseq[filter]` restriction was lifted for the remaining
seven groups:

| group | records | families | genera |
|---|---|---|---|
| anthozoa | 221 → 638 | 87 → 124 | 148 → 308 |
| porifera | 58 → 87 | 33 → 41 | 45 → 60 |
| mollusca | 850 → 1391 | 208 → 258 | 525 → 767 |
| arthropoda | 647 → 1078 | 169 → 231 | 358 → 561 |
| echinodermata | 135 → 263 | 56 → 80 | 104 → 171 |
| tunicata | 19 → 33 | 7 → 11 | 11 → 21 |
| annelida | 167 → 420 | 50 → 63 | 108 → 208 |

**Every rebuild is a strict superset**: checked by accession, no record that shipped in
the RefSeq-only build is missing from the new one, anthozoa's frozen 221 included. The
genus column is the one that matters — the point of lifting the filter is lineages that
had no representative at all, not more isolates of lineages that already did.

Anthozoa's 638 also absorbs the plain RefSeq refresh that used to be deferred here
(221 → 278 on RefSeq alone, from new entries and the widened rRNA-synonym matcher);
the new build holds all 278 of those plus 360 INSDC-only records.


The `.features.tsv` and `lineage` column were added to the then-current RefSeq-only
builds with `--refresh-derived` (below) rather than a rebuild, so that schema change
could not also change which records the databases held. `.fasta` and `.label.fasta` came
back byte-identical for all eight groups, and all 2,101 records of that generation
round-tripped through `refdb_record.py` to identical `coral_fix_bed.ref_features()` and
`reference_divergence_check.parse_reference()` output. That separation is why the record
counts above can be attributed to the widened search alone.

The completeness bar is per group and not negotiable upward for its own sake: coral
references must carry the features `CORAL_ANNOTATION_FIX` transfers (both rRNAs, a nad5
CDS) plus the 13-PCG cnidarian set, while ctenophore mitogenomes are genuinely reduced
(no atp6, no tRNAs, ~10 PCGs, rRNAs often unannotated), so the coral bar would reject
every valid ctenophore record.

## Why the RefSeq restriction was lifted

RefSeq is a curated subset, not a completeness bar — the records it omits are ordinary
INSDC submissions that pass exactly the same `record_is_complete()` check. Every group
was originally built with `AND refseq[filter]`, and for ctenophora that was the binding
constraint: 4 records out of 35 matching, with **no Platyctenida at all**. That is what
made `INV08_TJALFIELLA` look like an unfixable reference gap in the 20-sample invert
panel — it was seeded from three families in two other orders while **two *Tjalfiella*
mitogenomes**, its own genus, sat in GenBank behind the filter (PP327218, 11,397 bp,
11 CDS; PP331237, 11,020 bp, 11 CDS — both clear the group's own min-CDS bar of 10).

Ctenophora was lifted first, alone, because widening a group changes which record
`SELECT_REFERENCE_DB` picks and that was believed to re-pick the reference for corals
already submitted to ENA. **It does not.** `git ls-tree v2.0.0` has no `assets/refdb/`,
no `select_reference_db.py` and no `mito_origin_anchors.json`: the submitted corals were
assembled by the v2 path, which resolved a reference from the species *label* via
`findMitoReference`. No deposited mitogenome was ever built from these databases, and
every sample that has been through them is a test run. With the premise gone, the
remaining seven were lifted too.

What that bought, measured rather than assumed
([`bin/audit_reference_selection_diff.py`](../../bin/audit_reference_selection_diff.py),
50 assemblies from the invert panel, batch-20 and NOVA_260724_JP): **3 samples moved to
a closer reference, 0 moved further**, 13 changed record within the same taxonomic tier
and 34 kept the same record.

| sample | before | after |
|---|---|---|
| INV02_UMBELLULA | NC_044086.1 *Anthoptilum grandiflorum* (same order) | MK919668.1 *Umbellula huxleyi* (**congeneric**) |
| INV04_BOLOCERA | NC_066448.1 *Heteractis doreensis* (same order) | NC_022470.1 *Bolocera tuediae* (**congeneric**) |
| INV14_AMPHIOPHIURA | NC_085502.1 *Stegophiura sladeni* (same family) | LC698982.1 *Amphiophiura penichra* (**congeneric**) |

Two build-time dedup stages exist because of this, and are no-ops while the filter is on:

- **INSDC twins.** Without the filter a RefSeq record arrives alongside the
  submission it was derived from (`NC_038065` + `MG655622`). `drop_insdc_twins()`
  reads the `reference sequence is identical to X` comment and keeps the RefSeq copy.
  At scale this is load-bearing, not theoretical: 795 twins dropped in mollusca, 593 in
  arthropoda, 244 in anthozoa.
- **Per-organism cap** (`max_per_organism` in `GROUPS`, `--max-per-organism` to
  override). Set to 1 for the seven widened groups: a second isolate of a species adds
  no lineage the database did not already cover, and the tracked artifacts are plain git
  blobs rewritten wholesale on every rebuild. Ctenophora keeps 2, which is what its
  shipped 16-record build used — at that size a second isolate is worth its bulk. The
  cap removed 904 records in mollusca and 649 in arthropoda. Ranked RefSeq first, then
  longest, then most CDS, so a rebuild is reproducible.

Together those two stages are why lifting the filter costs far less than the raw search
counts suggest. Mollusca matches 3,169 records unfiltered but ships 1,391, and the eight
tracked databases went from 72 MiB to 133 MiB rather than the ~260 MiB a naive
record-count projection gives.

### Effect on the deposited-origin anchors

[`assets/taxonomy/mito_origin_anchors.json`](../taxonomy/mito_origin_anchors.json) is
derived from these databases, so a rebuild re-tallies it. Widening moved **no order off
its existing anchor**, and strengthened the evidence for the coral orders that carry the
already-deposited population — Scleractinia stays on trnM at n=58 → 181 (63.8% → 65.2%),
Malacalcyonacea on rrnL, Zoantharia and Scleralcyonacea on cox1. Two resolutions did
change, both for the better:

- **Actiniaria** now clears the bar in its own right (ND5, n=54, 79.6%) instead of
  falling back to the Anthozoa class aggregate and taking cox1.
- **Sabellida** drops below the bar and now resolves to its class (Polychaeta, cox1)
  rather than its old order-level trnH.

### Effect on coral annotation

`CORAL_ANNOTATION_FIX` is the one consumer that copies CONTENT out of the reference
rather than grading against it: `coral_fix_bed.py` BLAST-transfers the 16S rRNA and the
intron-split nad5 (and the cox1 exons when cox1 is intron-split too), so a changed
reference can change deposited sequence, not just a status line. "Same taxonomic tier"
does not establish that it does not.

Measured directly with
[`bin/audit_coral_annotation_diff.py`](../../bin/audit_coral_annotation_diff.py), which
re-annotates published corals from their own MITOS output (no MITOS2 re-run needed: the
reference enters only at `coral_fix_bed.py --ref-gb`, and `annotation/mitos_raw/result.bed`
is already published). 37 corals from batch-20 and NOVA_260724_JP, of which **11 receive a
different reference** under the widened database:

| samples | before | after |
|---|---|---|
| OG2327, OG2351, OG2354, OG2360 | *Hydnophora exesa* | *Coelastrea aspera* |
| OG2353, OG2357, OG2358, OG2359, OG2362, OG2367 | *Montipora efflorescens* | *Montipora mollis* |
| OG2352 | *Pavona decussata* | *Pavona clavus* |

**All 37 are byte-identical**, including those 11: same re-origined genome and same
per-gene translations, 13 PCGs, no internal stops, and the same `annotation_qc_gate.py`
verdict (35 PASS, 2 FIX) across every arm. The widened anthozoa database is inert for
coral annotation.

The audit runs three arms, because the published annotation was produced by older code on
this branch and a two-way diff could not tell a reference change from a code change:
`PUBLISHED` (on disk), `OLDREF` (today's code, old reference) and `NEWREF` (today's code,
new reference). `PUBLISHED` and `OLDREF` also agree for all 37, so coral annotation output
has not drifted either.

One implementation note worth keeping, because two obvious designs are wrong: the
cox1-rotated genome the MITOS BED refers to is **not** published (`annotation/<prefix>.fa`
is the re-origined genome, a different frame), so it is reconstructed and then verified
against the published annotation. It cannot be verified on cox1 position (`rotate_to_cox1.py`
does not always land cox1 at offset 0) nor on cox1 content (the fixer repairs cox1 -
OG2361's raw-BED cox1 is 873 bp against a published CO1 of 1572 bp). The check therefore
compares genes the fixer never rewrites, requiring several independent agreements.

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
