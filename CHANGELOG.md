# nf-core/oceangenomesmitogenomes: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

The invertebrate generalisation, on top of v2.0.0. Not released.

### `Added`

- The published origin of an invertebrate mitogenome is now measured per taxon instead
  of being tRNA-Met for everything.

  `mitos_to_emma.py` re-origined every invertebrate to trnM, under a docstring claiming
  that matched "how coral mitogenomes are deposited in NCBI". Measured against this
  repo's own curated RefSeq databases, trnM is the deposited origin for Porifera 0%,
  Annelida 0%, Ctenophora 0%, Echinodermata 0.7%, Mollusca 0.8% and Arthropoda 1.9%.
  Every one of the 12 annotated invertebrates in the 20-sample panel was published
  starting at `MT-TM`.

  The claim was true of exactly one lineage, the one it came from: **Scleractinia, at
  63.8% (37/58 records)**. That is also why the new table is keyed on taxonomic ORDER
  first, not class. Anthozoa as a class has no majority (rrnL 33.9%, cox1 28.5%, trnM
  17.6%) and that aggregate hides four orders with four different conventions:

  | order | n | anchor | share |
  |---|---|---|---|
  | Scleractinia | 58 | `TM` | 63.8% |
  | Malacalcyonacea | 75 | `RNR2` | 69.3% |
  | Zoantharia | 29 | `CO1` | 58.6% |
  | Scleralcyonacea | 31 | `CO1` | 51.6% |

  A class-level table would have rotated every already-submitted stony coral off its
  deposited origin and changed its ENA sequence checksum for no reason. Under the order
  level they are byte-identical: verified by re-running the adapter on OG2368's real
  MITOS output, which reproduces the published 17,842 bp sequence exactly.

  `bin/build_origin_anchor_table.py` generates `assets/taxonomy/mito_origin_anchors.json`
  by tallying which gene sits at position 1 across all 2,101 tracked reference records,
  resolving order/class via the taxdump. A taxon's own plurality wins when it clears both
  `min_fraction` (0.50) and `min_records` (20); otherwise it inherits order -> class ->
  group -> cox1. `InvertTaxonGroups.originAnchor(order, class)` reads it and never returns
  null: a wrong anchor rotates a circle, a missing one is a hard failure. The lookup is
  passed to MITOS2 and CORAL_ANNOTATION_FIX as `--origin-gene`, which is **required** in
  `mitos_to_emma.py` -- there is one rotation implementation with two call sites, and a
  default would let one be updated without the other, silently publishing FIX and PASS
  anthozoans of the same species on different origins.

  Notable non-cox1 results the phylum aggregates hid: Porifera `RNR2` 89.7% (Demospongiae
  92.7% at class level, Hexactinellida and Homoscleromorpha inheriting it), Echinoidea
  `TF` 75.6%, Cephalopoda `CO3` 73.8% (corroborated at order level by Sepiida 93.5% and
  Oegopsida 90.5%), Sabellida `TH` 68.0%.

  The pre-annotation cox1 rotation (`ROTATE_ORIGIN`) is unchanged and does a different
  job: it moves the *linearisation point* off an intron-split gene so MITOS annotates
  cleanly. The two cannot be merged, because tRNA-Met and rrnL are not findable by the
  tblastn protein search that step uses -- they need MITOS's annotation to exist first.

  **Operational constraint.** Changing an anchor changes the published sequence's
  checksum, so ENA treats a re-run sample as a different sequence even though the
  molecule is unchanged. Scleractinia keeps `TM`, so already-submitted stony corals are
  unaffected; the exposure is limited to submitted invertebrates whose anchor actually
  moves (octocorals, now `RNR2`, and sponges). Do not re-run those. If one genuinely
  needs resubmitting on its old origin, `originAnchorFor()` already prefers a per-sample
  `meta.origin_gene`, so pinning it needs only an `origin_gene` samplesheet column and no
  code change. There is deliberately no global "keep the old behaviour" flag: it would be
  the wrong default for every new sample and would persist silently.

- A GetOrganelle first pass that produces no contig at all can now be reseeded.

  `INV10_CHITON` and `INV12_ANOMURA` both recruited reads and built a graph, then logged
  `Slimming ... finished with no target organelle contigs found!` / `No sequence hit our
  LabelDatabase!` and published a zero-byte assembly. That is a *labelling* miss, not a
  recruitment miss: the stock `animal_mt` LabelDatabase cannot label a divergent
  invertebrate contig. Representation was never the problem -- the curated databases hold
  20 Polyplacophora records (6 of them Chitonidae) and 20 Anomura across 9 families.

  Two changes. First, the first pass is now labelled with the sample's curated group gene
  database (`--genes`), resolved from `meta.class` alone via
  `InvertTaxonGroups.seedDbGroup()` with no assembly involved, which is what makes it
  available on the first pass at all. The stock seed (`-F animal_mt`) still recruits the
  reads, so recruitment is unchanged by construction. A vertebrate, or a class with no
  curated database, passes `[]` and emits a byte-identical command line. Revert with
  `--getorganelle_firstpass_group_genes false`.

  Second, an empty first pass is reseeded from the WHOLE group database. Until now the
  reseed's seed was chosen by BLASTing the first-pass assembly against the group, so a
  total first-pass failure -- the case that most needs a better seed -- could never
  select one: `select_reference_db.py` reported `NONE  empty assembly`, readiness came
  back false, and the sample was routed straight back to the empty assembly it already
  had. Coarse whole-group seeding is what the top-n subset exists to avoid, so it fires
  only on a zero-byte assembly, where `preferReseed()` cannot make the outcome worse.
  Revert with `--getorganelle_empty_first_pass_rescue false`.

- Samples that leave during ENA preparation are now reported as held instead of vanishing.

  Four QC-passing invertebrates (INV03_ZOANTHUS, INV11_EUCRATE, INV13_SEMIBALANUS,
  INV15_ACANTHASTER) appeared in neither `held_samples.tsv` nor the submission-ready set,
  and `run_completeness.txt` could only report them as unaccounted. `BUILD_SOURCE_MODIFIERS`
  had exited 1; it has no explicit `errorStrategy`, so it inherited the global `ignore`,
  emitted no `.src`, and the inner join at `ch_processed_files` dropped the sample silently.

  `ENA_SUBMISSION_PREP` now emits a `PRE_SUBMISSION` held fragment for the remainder of its
  own input that never reached the table2asn verdict, naming the stage that produced no
  output. A remainder join rather than an `errorStrategy` change on that one process:
  catching only `BUILD_SOURCE_MODIFIERS` would leave every sibling step with the same hole.

- Invertebrate reseeds now narrow the group database down to the records the sample's
  own assembly matches, instead of seeding from the whole thing.

  Resolving a seed is two narrowing stages. `InvertTaxonGroups.seedDbGroup()` (above) added
  the first, by taxonomy: class -> one of the eight group databases. The second, by
  sequence, existed only for vertebrates, where `findMitoReference` resolves a single
  relative and `REFERENCE_RANK` re-picks it from the reads. Invertebrates stopped after
  stage 1 and handed GetOrganelle the entire group -- 221 anthozoan genomes across 87
  families, and their 3,776 gene sequences, as `-s` and `--genes`. That recruits reads from
  across the phylum: INV04_BOLOCERA's reseed converged on 46,142 accepted words against its
  first pass's 16,192, average base-coverage fell from 26.0 to 11.7, and a 2-scaffold first
  pass came back as 12 scaffolds. Anthozoa is the *mildest* case that can break this way --
  mollusca is 850 genomes and 31,718 label sequences, arthropoda 647 and 23,991.

  `SELECT_REFERENCE_DB` (was `SELECT_CORAL_REFERENCE`, and now group-agnostic) BLASTs the
  first-pass assembly against the group and hands the reseed only the top
  `params.reseed_seed_top_n` records, with just their genes: 88 label sequences rather than
  3,776 for INV04_BOLOCERA. Top-n rather than the single best, because from a small
  fragmented first pass the pick is reliable at order/subclass level but not at species --
  INV01_ACANTHOGORGIA's own family is in the database and does not win from 1,194 bp, while
  the top 5 are all Octocorallia. A sample nothing aligns to is **not** reseeded: that is
  the honest signal that the group does not represent its lineage, and falling back to the
  whole group is the failure being removed.

- Invertebrates now carry a real reference through `GETORGANELLE_CHECK`.

  Every invertebrate previously reached the check with the `NO_REFERENCE.gb` placeholder, so
  its evidence row recorded `note=no_reference` with `NA` reference coverage and length
  ratio -- INV04_BOLOCERA's 12-scaffold reseed produced no length or coverage signal at all.
  The seed-mode selector emits the single best record alongside the top-n seed, off the same
  ranking and at no extra BLAST, and that becomes the reference. It also lets
  `REFERENCE_RELEVANCE` reach the invertebrate branch it already had (`min_pid` 88 vs 82) and
  which no sample could previously satisfy.

- Invertebrates whose first GetOrganelle pass produced no contig are now reseeded from a
  taxonomy-bounded, read-ranked panel instead of the whole group database
  (`SELECT_FALLBACK_SEED`, `bin/select_fallback_seed.py`).

  An empty first pass has no sequence for `SELECT_REFERENCE_DB` to rank a group database
  against, so the rescue handed GetOrganelle the entire group. That is defensible for a
  4-record database and indefensible for a 647-record one: `INV12_ANOMURA` was seeded from
  all of Arthropoda while **20 Anomura records sat in the same manifest**. Whole-phylum
  seeding recruits conserved and off-target reads, which is exactly what stage 2 exists to
  prevent.

  The new selector walks the manifest lineage nominal → family → order → class and stops at
  the most specific rank that returns records, then ranks that shortlist by mapping a fixed
  read subsample. Against the tracked databases: Anomura resolves 20 of 647 at tier `order`,
  Chitonida 17 of 850, and `Tjalfiella sp.` resolves 2 records at tier `nominal`. When no
  rank matches it still falls back to the group, but names every rank that missed
  (`reference_gap=nominal:…,family:…,order:…`) so a genuine reference gap is auditable
  rather than indistinguishable from an ordinary selection. `INV05_FARREA` is the live
  example: *Farrea* has no complete mitogenome in GenBank at all, so no amount of ranking
  can help it and the run should say so.

  Bounded rather than best-only for the same reason the seed selector is top-n: from no
  assembly at all, the pick is trustworthy at order level and not at species. `balanced_cap`
  spreads a large tier across families so the panel stays broad within its bound.

- The reference-database builder takes `refseq_only` per group, and ctenophora is the first
  group to lift it: **4 records → 16, 3 families → 7**.

  `SEARCH_TEMPLATE` applied `AND refseq[filter]` to every group. RefSeq is a curated subset,
  not a completeness bar — the records it excludes are ordinary INSDC submissions that pass
  the identical `record_is_complete()` check. For seven groups the restriction is a harmless
  de-duplicator. For ctenophora it was the binding constraint and it hid most of the phylum:
  4 records out of 35 matching, **and no Platyctenida whatsoever**.

  That is what made `INV08_TJALFIELLA` look like an unfixable reference gap. It was not. Two
  *Tjalfiella* mitogenomes — its own genus — were sitting in GenBank behind the filter, and
  both clear the group's own bar of 10 CDS:

  | accession | bp | CDS | rRNA | organism |
  |---|---|---|---|---|
  | PP327218 | 11,397 | 11 | 2 | *Tjalfiella* sp. |
  | PP331237 | 11,020 | 11 | 2 | *Tjalfiella* sp. |

  The rebuild adds Tjalfiellidae, Lyroctenidae, Benthoplanidae and Euplokamidae. `min_cds`
  is untouched at 10 — the platyctenids clear it, and the two 9-CDS *Pleurobrachia* records
  are still correctly dropped.

  It was lifted for ctenophora alone at first, on the belief that widening the other seven
  would re-pick the `SELECT_REFERENCE_DB` reference for samples already submitted. That
  belief was wrong and the remaining seven were lifted immediately after — see the next
  entry. Two dedup stages exist only because of this and are no-ops while the filter is on: `drop_insdc_twins()` removes the INSDC
  submission a RefSeq record was derived from (`NC_038065` + `MG655622`), and
  `--max-per-organism` (default 2) stops one heavily-resequenced species filling the panel —
  the widened ctenophore search returns **nine** *Vallicula multiformis* isolates, all in one
  family, where neither a top-n panel nor family balancing could dilute them.

  `assets/taxonomy/mito_origin_anchors.json` is regenerated with it. Only ctenophore tallies
  move and every one stays below `min_records=20`, so no anchor is gained or lost and no
  submitted sequence is re-origined. The staleness warning now fires for a single-group
  rebuild too, not only `--all`, which is how this was nearly missed.

- Every remaining group now searches all of INSDC as well, at one record per organism.
  **2,101 records → 3,926 across the eight, 617 families → 815, and 1,308 genera → 2,105.**

  | group | records | families | genera |
  |---|---|---|---|
  | anthozoa | 221 → 638 | 87 → 124 | 148 → 308 |
  | mollusca | 850 → 1391 | 208 → 258 | 525 → 767 |
  | arthropoda | 647 → 1078 | 169 → 231 | 358 → 561 |
  | annelida | 167 → 420 | 50 → 63 | 108 → 208 |
  | echinodermata | 135 → 263 | 56 → 80 | 104 → 171 |
  | porifera | 58 → 87 | 33 → 41 | 45 → 60 |
  | tunicata | 19 → 33 | 7 → 11 | 11 → 21 |

  The reason the other seven were held back does not exist. `git ls-tree v2.0.0` has no
  `assets/refdb/`, no `select_reference_db.py` and no `mito_origin_anchors.json`: the
  corals already in ENA were assembled by the v2 path, which resolved a reference from the
  species *label* through `findMitoReference`. No deposited mitogenome was ever built from
  these databases, and every sample that has been through them is a test run.

  Measured rather than assumed, with the new `bin/audit_reference_selection_diff.py` over
  50 published assemblies from the invert panel, batch-20 and NOVA_260724_JP: **3 samples
  moved to a closer reference and 0 moved further**, 13 changed record inside the same
  taxonomic tier and 34 kept the same record.

  | sample | before | after |
  |---|---|---|
  | INV02_UMBELLULA | NC_044086.1 *Anthoptilum grandiflorum* (same order) | MK919668.1 *Umbellula huxleyi* (**congeneric**) |
  | INV04_BOLOCERA | NC_066448.1 *Heteractis doreensis* (same order) | NC_022470.1 *Bolocera tuediae* (**congeneric**) |
  | INV14_AMPHIOPHIURA | NC_085502.1 *Stegophiura sladeni* (same family) | LC698982.1 *Amphiophiura penichra* (**congeneric**) |

  Every rebuild is a **strict superset**, checked by accession: no record that shipped in
  the RefSeq-only build is absent from the new one, anthozoa's deliberately frozen 221
  included. Anthozoa's 638 also absorbs the plain RefSeq refresh that was being deferred
  (221 → 278 on RefSeq alone), so the freeze is retired rather than merely overtaken.

  `max_per_organism` joins `refseq_only` in `GROUPS` and is 1 for the seven, so a rebuild
  is reproducible from `--group X` with no flag to remember. That and `drop_insdc_twins()`
  are what keep the cost down: mollusca matches 3,169 records unfiltered but ships 1,391
  (795 INSDC twins and 904 over-cap records removed), and the tracked databases grow
  72 MiB → 133 MiB rather than the ~260 MiB the raw counts imply. Ctenophora keeps its cap
  of 2 and its 16 records are untouched. `--retmax` had to rise from 2000 to 6000, since
  `fetch_records()` aborts rather than silently shipping an arbitrary slice of a group.

  `assets/taxonomy/mito_origin_anchors.json` is regenerated with the rebuild and **no order
  moved off its anchor**. The coral orders that carry the deposited population only gain
  evidence: Scleractinia stays on `TM` at n=58 → 181 (63.8% → 65.2%), Malacalcyonacea on
  `RNR2`, Zoantharia and Scleralcyonacea on `CO1`. Two resolutions do change, both toward the
  measured answer: Actiniaria now clears the bar in its own right (`ND5`, n=54, 79.6%)
  instead of taking `CO1` from the Anthozoa class aggregate, and Sabellida falls below it
  and resolves to its class (Polychaeta, `CO1`) rather than its old order-level `TH`.

  Timing was checked before assuming a cost: `SELECT_REFERENCE_DB` runs `blastn -subject`
  with no formatted database, and a mollusc selection takes 3 s against both the 13 MiB and
  the 22 MiB subject, so no `makeblastdb` step or resource-label change is warranted.

### `Fixed`

- `DOWNLOAD_BLAST_DB` uses the shared taxdb cache again, instead of re-downloading
  62 MB over unguarded FTP on every single run.

  Three faults had stacked up. The `storeDir` cache at `--blast_db_dir` was
  permanently unsatisfiable: `taxonomy4blast.sqlite3` was added as a required output
  without seeding the existing cache, so Nextflow logged `unable to find ...
  taxonomy4blast.sqlite3` and re-submitted the task every run. Nothing consumes that
  output, and blastn only ever receives `taxdb.btd`/`taxdb.bti`, so it is now
  `optional: true` and a two-file cache is a hit.

  The in-script guard that was meant to short-circuit the download had been dead since
  the cache base moved from `${launchDir}/blast_dbs` to `params.blast_db_dir`: the
  literal `/blast_dbs` suffix stayed behind, so it tested
  `<blast_db_dir>/blast_dbs/taxdb.btd`, a path that never existed. The suffix is gone,
  and on a hit the guard now *copies* the cached files into the work dir rather than
  echoing "skipping" -- a bare skip would leave the declared `taxdb*` outputs
  unresolvable, because `storeDir` runs the task in a clean work dir.

  Finally the download itself was unverified. A dropped FTP control connection resumed
  via `REST`, reported all 65,261,809 bytes saved, and produced a corrupt archive;
  `tar` exited 2, which is deliberately outside the transient set in `base.config`, so
  the whole cohort terminated. The fetch is now https (no control/data channel to
  drop), retries, and checks NCBI's `taxdb.tar.gz.md5` before unpacking, deleting the
  partial file on mismatch so a retry cannot resume onto known-bad bytes.

  `MITOGENOME_ANNOTATION` now also errors up front when `--blast_db_dir` is unset,
  mirroring the existing `--taxonkit_db_dir` assert, rather than silently disabling the
  cache via `storeDir null`.

- BLAST summaries reach MultiQC again, instead of 48 copies of a directory name.

  `BLAST_BLASTN` emits `summary` as a plain path and `tool_params` as a `[meta, path]`
  tuple, on adjacent lines of the same module. The MultiQC collection applied the tuple
  idiom to the plain one:

  ```groovy
  ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect { it[1] })
  ```

  Indexing a `Path` does not fail — it returns the name element at that index. For a summary
  under `/scratch/pawsey1348/…` element 1 is the literal string `pawsey1348`, so Nextflow
  tried to hash and stage a directory component as a MultiQC input, once per summary:

  ```
  WARN: [HashBuilder] Unable to get file attributes file: /<outdir>/pawsey1348
  ```

  Every real `filtered_summary.*.txt` was dropped in the process, and the run reported only
  a warning. Now mixed directly. `tests/multiqc_inputs/` runs both idioms over one channel,
  so the fix is pinned by what it produces rather than by the shape of the source line.

- `-profile test` runs. It had never run since the repository was scaffolded.

  `conf/test.config` still carried the nf-core template placeholders, `// TODO nf-core:
  Specify the paths to your test data` included: `input` pointed at an nf-core *viralrecon*
  amplicon samplesheet and `genome` at the yeast assembly `R64-1-1`. That samplesheet has no
  `sequencing_type` column, which `assets/schema_input.json` marks required, so every
  invocation died in parameter validation before a single task was scheduled:

  ```
  * --input (.../viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv):
      Entry 1: Missing required field(s): sequencing_type
  ```

  The profile now uses `tests/test_data/samplesheet_test.csv`, two invertebrate HiFi rows
  carrying the full taxonomy columns, and sets `organelle_type = 'animal_mt'` (null by
  default, and `GETORGANELLE_CONFIG` takes it as a plain `val`, so the run aborted at
  channel construction once the samplesheet was fixed). The annotation, ENA prep and upload
  stages are skipped: each hard-requires an external resource this repository cannot ship
  (`--taxonkit_db_dir`, `--nt_blast_db`, `--template_sbt`, Postgres credentials). What
  remains is a stub smoke test of samplesheet parsing and the assembly wiring.

  The reads are bundled, but nf-schema resolves a relative `fastq_1` cell against
  `launchDir`, not against the samplesheet, so a committed relative path would only work
  when launched from the clone root. The tracked sheet keeps the `__PROJECT_DIR__`
  placeholder that the nf-tests under `tests/samplesheet_meta/` already substitute, and the
  profile absolutises it into a run-local temp copy. Verified from both the clone root and
  an unrelated working directory.

- `reseed_seed_top_n` and `reseed_length_tolerance` are declared in `nextflow_schema.json`.

  Both had `nextflow.config` defaults but no schema entry, so every run opened by reporting
  its own defaults as invalid. Declared with bounds (`top_n` 1–20, `tolerance` ≥ 1.0)
  alongside the two new `reseed_fallback_*` parameters.

### `Changed`

- `GETORGANELLE_FROMREADS`/`GETORGANELLE_RESEED`/`OATK` request 64 GB on their first
  attempt instead of 100 GB (`params.hinted_memory_base_gb`).

  Setonix bills `max(cpus/128, memory/230) * 128` core-equivalents, so the memory request
  is what is charged: 100 GB costs 55.7 cores against the 16 these processes actually ask
  for. That is why `GETORGANELLE_FROMREADS` alone was 548.63 of the 20-sample invert
  panel's 651.26 SU (84%), at 3.14% median CPU and 13.02% median memory efficiency. 64 GB
  charges 35.6 cores, a 36% cut, and is the lowest round tier above the 51.9 GB peak RSS
  observed across those 19 assemblies (min 6.9 GB, median ~13 GB); the retry ladder
  (128, 192 GB) covers outliers. `cpus` is deliberately unchanged -- 3.14% is real, but
  cpus and walltime interact with GetOrganelle's own threading and need a benchmark first.

  Resource directives are not hashed, so this does not invalidate `-resume`. It does not
  reach a sample that already has a per-sample hint in `pipeline_info/memory_hints.json`,
  though: hints are only recorded on a retry, so every existing one is the old 200 GB tier
  and will replay verbatim. **Operational step:** drop the `GETORGANELLE_FROMREADS` and
  `GETORGANELLE_RESEED` blocks from that file in each production `outdir` so the new ladder
  re-derives them. Not automated on purpose -- that file is evidence of what actually
  succeeded.

- The reseed replaces the first-pass assembly only when it is actually better.

  The resolution was `rs_fasta.size() > 0` under a comment claiming it kept "the better of"
  the two, so any non-empty reseed won -- which is how INV04_BOLOCERA published a 12-scaffold
  reseed over its own 2-scaffold first pass. `preferReseed()` now ranks them: an empty reseed
  still never wins; circularity decides first, on the same log evidence `needsReseed` uses;
  then a length difference beyond `params.reseed_length_tolerance` (1.5x), because a reseed
  that recovered several times more sequence in more pieces is not a regression; then
  contiguity. The loser is unaffected in the database -- `selectProvenanceVariants()` already
  emits every superseded attempt as its own row.

- `assembly_length` in `*.getorg_check.tsv` is the whole assembly, not its first record.

  `check_getorganelle.py` computed the total and then reported `len(first_record)`, so
  INV04_BOLOCERA's 12-scaffold, 14,745 bp reseed was recorded as 268 bp. The first record is
  still what the reference-coverage, length-ratio and tandem-repeat tests run on -- they
  describe one molecule and must not be fed a concatenation -- and is now reported as its own
  `first_record_length` column, appended at the end. `bin/check_circularity.py` (OATK) keeps
  first-record semantics deliberately: there a second record is an anomaly it warns about,
  not the expected shape.

- The per-group `.gb` files are no longer tracked, for any group.

  They were 84 MB raw / 26 MB compressed across the eight, plain git blobs with no LFS, and
  every rebuild rewrote all of them into history -- which is why only anthozoa's was tracked,
  and therefore why sequence-based reference selection was confined to corals. Everything the
  pipeline reads out of a reference GenBank is now in four small tracked files: sequence and
  length from `.fasta`, organism and taxonomy from a new `lineage` column on
  `.manifest.tsv`, and feature types with their exon coordinates from a new
  `.features.tsv` (3.6 MB raw / 0.7 MB compressed for all eight). `bin/refdb_record.py`
  rebuilds an equivalent record on demand, so no downstream consumer changed: they still take
  a `--ref-gb`/`--reference-gb` path, which is also what keeps them working for vertebrates,
  whose reference really is a GenBank downloaded by `findMitoReference`.

  `.features.tsv` cannot be replaced by `.label.fasta`: the label database stores each
  feature's *spliced* sequence (`feat.extract`), while `coral_fix_bed.py` needs one entry per
  exon of the group-I-intron-split nad5. `build_invert_reference_db.py --refresh-derived`
  regenerates the derived files from an existing `.gb` without contacting NCBI, so a schema
  change cannot smuggle in a content change -- used here, and `.fasta`/`.label.fasta` came
  back byte-identical for all eight groups, with all 2,101 records round-tripping to
  identical `ref_features()` and `parse_reference()` output.

### `Added`

- Invertebrate GetOrganelle reseeds now seed from a curated database for the sample's own
  phylum instead of the coral one.

  A failed first pass was reseeded from `coral_mito_refdb.*` for *every* invertebrate, and
  that database is Anthozoa-only (`txid6101`). A mollusc or a sea star was therefore re-run
  against a seed far too divergent to assemble from: a guaranteed second failure that still
  cost the full GetOrganelle walltime. Seven new databases are built alongside it, one per
  group we run -- porifera, mollusca, arthropoda, echinodermata, ctenophora, tunicata,
  annelida -- and `InvertTaxonGroups.seedDbGroup()` resolves each sample's class to one of
  them. A class no database covers now returns null and is **not** reseeded at all: it keeps
  its first-pass assembly and is named in the log, because not reseeding is strictly better
  than reseeding from the wrong phylum. This is deliberately the opposite of
  `cox1PanelGroup()`'s catch-all -- a mismatched rotation panel is a safe no-op, a
  mismatched seed is a wasted run.

  `ARTHROPODA_CLASSES` gains Thecostraca, Copepoda, Ichthyostraca and Mystacocarida (the
  classes the old Maxillopoda was split into, all of which NCBI returns), and each phylum set
  gains its own phylum name for samples whose taxonomy only resolved to phylum rank. Without
  them a barnacle got no seed database and the *coral* cox1 rotation panel. New
  `tests/invert_seed_db` resolves every class we run and asserts the files are on disk, since
  the group name is also the directory name -- the drift that hid the barnacle gap.

  `bin/build_coral_reference_db.py` becomes `bin/build_invert_reference_db.py`, taking
  `--group`/`--all` over a registry of Entrez organism expressions and per-group completeness
  bars. The bar cannot be one number: coral references must carry what `CORAL_ANNOTATION_FIX`
  transfers (both rRNAs, a nad5 CDS) plus the 13-PCG cnidarian set, while ctenophore
  mitogenomes are genuinely reduced (no atp6, no tRNAs, ~10 PCGs), so the coral bar rejects
  every valid ctenophore record. The rRNA matcher also learns the `l-rRNA`/`rnl`/`MT-RNR2`
  spellings other phyla use, where the coral-era `16S|RRNL|LARGE` test silently dropped
  correctly annotated records; Arthropoda subtracts Hexapoda, Arachnida and Myriapoda, which
  are >95% of arthropod RefSeq mitogenomes, none of them anything OceanOmics sequences, and
  enough to truncate the search at `--retmax`. The builder now aborts on such a truncation
  rather than shipping an arbitrary slice of a group.

  The anthozoa database is kept at its existing 221-record build. Rebuilding it today yields
  278 -- a strict superset, no record lost -- but those extra records also change which
  reference `SELECT_CORAL_REFERENCE` picks for every coral, so that refresh belongs in its
  own change.

- `assets/` is reorganised from 35 flat files into `refdb/<group>/`, `panels/`, `taxonomy/`
  and `placeholders/`, with nf-core boilerplate left at the root. Only the files the pipeline
  reads are tracked: the per-group `.gb` files (85 MB, read by nothing today) and the
  `.nucl.*` BLAST databases are gitignored and rebuilt on demand. `assets/refdb/README.md`
  records each database's taxon, completeness bar, record count and rebuild command -- the
  provenance that previously existed only in the builder's docstring.

- Invertebrate annotation now generalises beyond Cnidaria to Mollusca, Echinodermata,
  Arthropoda (Crustacea, Pycnogonida) and Porifera, instead of running every invertebrate
  through coral-tuned machinery. `mitoGeneticCode()` gains an explicit allow-list of the
  invertebrate classes confirmed to use the standard Invertebrate code (NCBI 5, previously
  everything invertebrate took the Coelenterate code 4) and a Porifera Coelenterate-code (4)
  entry alongside Cnidaria; there is no catch-all code-5 default -- an invertebrate class in
  none of the code-4, code-9 or code-5 groups aborts the run so its code is resolved before
  annotation rather than guessed, since not every invertebrate is code 5 (Bivalvia among
  others) and a wrong table mistranslates every CDS and fails table2asn terminally. The
  per-sample `genetic_code` samplesheet column overrides the map;
  `ROTATE_ORIGIN` re-origins each sample against its own phylum-appropriate curated cox1 panel
  (`assets/panels/cox1/{mollusca,arthropoda,echinodermata}.faa`, new) instead of the coral-only
  `panels/cox1/anthozoa.faa` for every invertebrate; `REFERENCE_RELEVANCE`'s 88.0 identity floor is now
  Cnidaria-only, with everything else defaulting to the vertebrate 82.0 as an untuned starting
  point. `ANNOTATION_QC_GATE`/`CORAL_ANNOTATION_FIX` stay Cnidaria-only rather than guessing an
  equivalent for other phyla -- the nad5 group I intron they repair is a Hexacorallia-specific
  trait, absent in Octocorallia and not the same failure mode Porifera's own (cox1, lineage-
  specific) group I introns present; non-Cnidarian invertebrate MITOS2 output now merges straight
  through untouched. New `InvertTaxonGroups` (`lib/`) centralises the class groupings both the
  samplesheet and annotation subworkflow read.

- The generated samplesheet now carries a resolved `genetic_code` per row, and the class ->
  code map lives in one place.

  `bin/create_samplesheet.py` emitted the `genetic_code` column but always wrote it blank,
  leaving `mitoGeneticCode()` to derive the code from the class. That was harmless while
  `mitoGeneticCode()` had a catch-all code-5 default for invertebrates; once it aborted on an
  unmapped class instead, a blank column meant the run died at samplesheet parse for any class
  the map did not know. It knew 25 of the 76 classes `INVERT_CLASSES` marks as invertebrate:
  Bivalvia, Cephalopoda, Polychaeta, Ascidiacea and 47 others would have stopped the run.

  New `assets/taxonomy/mito_genetic_codes.json` is the single source of truth, read by both
  `bin/create_samplesheet.py` and `lib/InvertTaxonGroups.groovy` -- the same fix
  `bin/mito_gene_order.py` applied to `REF_GENES` after four copies drifted. Every entry
  records its `basis` and whether the code is NCBI-documented (`ncbi`) or the conventional
  choice for the group (`convention`), and an `ambiguous` block documents the classes that are
  deliberately left to abort: Pterobranchia and Hemichordata (Rhabdopleuridae is code 24,
  Cephalodiscidae 33, unresolvable at class rank), Trematoda/Cestoda/Monogenea/Platyhelminthes
  (NCBI table 21 is trematode-specific and 14 is the alternative flatworm code, so table 9's
  stated Rhabditophora scope does not cover them -- these were previously mapped to 9) and
  Placozoa. Ctenophora joins Cnidaria on code 4, which is NCBI table 4's stated metazoan scope
  verbatim; tunicates take the Ascidian code 13.

  The generator resolves each row against that map and writes the code, so the sheet shows
  which table a sample will be annotated under and can be corrected before launching.
  Unresolved invertebrate rows are still written, and are now named on stderr at generation
  time (`report_unresolved_genetic_code()`, mirroring the existing `report_unresolved()`),
  with the ambiguity reason where there is one, instead of surfacing as an abort mid-run.
  `taxonomy_resolution.tsv` records the code alongside the taxonomy provenance.

- The genetic code now comes from NCBI's own per-taxon assignment, not from a class map.

  `nodes.dmp` field 8 holds the mitochondrial genetic code NCBI assigns each taxon, and it is
  what ENA and table2asn validate a submission against. `TaxdumpLineage` now parses it and
  `resolve_genetic_code()` prefers it whenever the taxdump resolved the sample, with
  `assets/taxonomy/mito_genetic_codes.json` demoted to the fallback for a class that came from the
  species table with no lineage behind it.

  Per-taxon beats per-class because the code varies *below* class rank: `Cephalodiscidae` is
  code 33 while its parent class `Pterobranchia` is 5, so no table keyed on class can hold
  both. A `Cephalodiscus` sample now gets 33 and a `Rhabdopleura` 5, from the same class.

  Checking the hand-written map against that field found two errors in it: `Appendicularia`
  (larvaceans, e.g. `Oikopleura`) and the `Tunicata` subphylum node are **code 5, not 13**.
  NCBI's prose page scopes table 13 as "Urochordata (tunicates)", but its taxonomy applies 13
  only to Ascidiacea and Thaliacea. The remaining 70 entries agree with NCBI exactly, and a
  new test asserts that agreement so the map cannot drift from it again.

  With per-taxon codes available, the `ambiguous` block shrinks to the pterobranchs, the only
  group whose code genuinely varies below class. `Trematoda`, `Cestoda`, `Monogenea`,
  `Platyhelminthes` (all 9) and `Placozoa` (4) are mapped rather than left to abort -- NCBI is
  unambiguous on them, and the earlier caution was based on the prose page's table 14/21
  scoping rather than on what NCBI assigns.

- Samples identified only to phylum are now called out at samplesheet-generation time.
  `report_phylum_only()` names them on stderr with the phylum, alongside the existing
  unresolved-taxonomy and unresolved-code reports, so a phylum sitting in a column called
  `class` is stated rather than left to be discovered by whoever reads the sheet next.

- Invertebrate status is now decided by NCBI ancestry, not by a hand-maintained class list.

  `INVERT_CLASSES` is an allow-list, and anything missing from it read as a vertebrate:
  genetic code 2, EMMA instead of MITOS2 and the curated fish BLAST database, with nothing
  raised to notice. Checking the first invertebrate run against the taxdump found stalked
  barnacles going through as fish -- they resolve to class `Thecostraca`, which was absent --
  and a sweep of every class-rank node under Metazoa but outside Vertebrata put **41 of NCBI's
  92 invertebrate classes** outside the list, Insecta, Arachnida, Copepoda and Merostomata
  among them. `Maxillopoda` and `Hexanauplia`, which were in the list, are retired names the
  taxdump no longer carries at all.

  `TaxdumpLineage` now reports `is_animal` / `is_vertebrate` from the lineage walk, and
  `is_invertebrate()` uses them whenever the taxdump resolved the sample, falling back to
  `INVERT_CLASSES` only for a class that came from the species table with no lineage behind it.
  A class nobody has added can no longer misroute a sample. The ancestry anchors are picked by
  lineage rather than by first match, because `Vertebrata` is itself a homonym (a red algal
  genus) and taking the wrong one would classify every vertebrate as an invertebrate.

- A cross-kingdom homonym now resolves to the animal rather than to nothing. Dropping the name
  outright was too blunt for a pipeline that sequences animals and never plants: it cost the
  sponge genus `Acanthella` and the barnacle genus `Calantica` their whole lineage, both being
  homonyms of flowering-plant genera. When exactly one candidate is an animal it is the answer;
  two animal candidates (the phylum `Ctenophora` and the crane-fly genus `Ctenophora`) stay
  genuinely ambiguous and still resolve to nothing.

- A sample identified no further than its phylum now resolves. `phylum` joins `INDEXED_RANKS`
  and `WANTED_RANKS`, and a lineage with no class rank falls back to the phylum name as the
  sample's `class`. A sponge recorded only as `Porifera` pins no class, but the phylum still
  selects the genetic code (4), the reduced-tRNA expectation and the cox1 rotation panel --
  `INVERT_CLASSES` and `mito_genetic_codes.json` already carry the phylum names for exactly
  this case. These rows are marked `taxdump-phylum` in `taxonomy_resolution.tsv` rather than
  being hidden under `taxdump`.

- Barnacles are no longer annotated as vertebrates: `Thecostraca` and `Copepoda` (plus
  `Ichthyostraca` and `Mystacocarida`) join the crustacean classes, all on code 5.

- `Craniata` is no longer treated as an invertebrate class. It is a brachiopod class, but it is
  also the vertebrate clade name, so a fish whose class resolved to `Craniata` was marked
  `invertebrates=true` and routed to the invertebrate BLAST DB and MITOS2. Craniate brachiopods
  use `Rhynchonellata`/`Lingulata`/`Brachiopoda`, or the `genetic_code` column.

- `MITOS2`'s stub block no longer falls back to genetic code 5 when `meta.genetic_code` is
  unset, matching the script block. Both now rely on the upstream abort and the
  `SUPPORTED_GENETIC_CODES` assertion rather than guessing a table.

### `Fixed`

- ENA submission prep no longer requires the database *writes*, only a database. It now
  runs whenever `--sql_config` is supplied, including under `--skip_upload_results`.

  Narrowing `--skip_upload_results` to mean "skip the PostgreSQL writes" (below) freed local
  QC but deliberately left `ENA_SUBMISSION_PREP` inside the write guard, on the grounds that
  it needs the database. It does -- but only to read: `BUILD_SOURCE_MODIFIERS` and
  `PREPARE_ENA_METADATA` want collection date, country, coordinates and coverage. Supplying
  `--sql_config` was therefore not enough to make it run, which is not what the flag says.

  One thing actually coupled it to the writes. `PREPARE_ENA_METADATA` read `mean_depth` back
  out of the `mitogenome_data` row that this same run had just written, so submission prep
  had to be ordered behind the committed row with `gateOnAssemblyReceipt`, and with the
  writes off there were no receipts for the gate to join against. That round trip bought
  nothing: `MITOGENOME_COVERAGE` measures the depth, `push_mtdna_assm_results.py` writes that
  number to the row, and this script read the same number back. It also had a quieter cost --
  on a re-run that skipped the write it returned the *previous* assembly's depth for a
  molecule that had just been reassembled.

  Coverage now travels to submission prep as the depth TSV itself, carried on the QC tuple
  and attached with the same remainder-join-plus-placeholder pattern as the assembly upload
  rows. `bin/depth_tsv.py` reads it without pandas (the metadata container has psycopg2 and
  not much else); the database `SELECT` survives as the fallback for runs that measured no
  depth, and `mean_depth_source` in the package metadata records which one answered.
  `gateOnAssemblyReceipt` is untouched and still guards the four pushers in `UPLOAD_RESULTS`.

  `--skip_ena_submission_prep` is the new off switch for the stage, since
  `--skip_upload_results` no longer doubles as one. Prep also now requires `--ena_study`,
  which has no default: it hard-errors without one, and that error was unreachable while
  prep sat behind the upload guard. Rather than abort every upload-only run that never
  intended to submit, the study gates prep -- absent, it is skipped with a warning; present,
  prep keeps its strict PRJEB accession check. `nextflow_run_external.sh` passes
  `--skip_upload_results true`, so it will now run submission prep whenever it is also given
  a `--sql_config`.

  One related silent failure closed on the way: `build_source_modifiers.py` SELECTs from
  `sample` but filters on `mitogenome_data` (`WHERE m.og_id`/`m.tech`), and an empty result
  wrote two empty CSVs and exited 0. That was survivable only while the receipt gate
  guaranteed the parent row existed. Prep now runs in modes where a sample genuinely has no
  row, and an empty `.src` would reach a submitter looking exactly like a specimen with
  nothing recorded about it, so it fails loudly instead.

- QC no longer requires a database. `--skip_upload_results` now means "skip the PostgreSQL
  writes" and nothing else.

  Everything that answers "is this mitogenome any good?" lived inside
  `if (!params.skip_upload_results && params.sql_config)`. A run with `--skip_upload_results`,
  or with no `--sql_config` at all -- which is how the invertebrate work is developed --
  therefore produced no gene counts, no missing-gene set, no completeness verdict, no
  held-samples report, and a `mitogenome_assembly_summary_mqc.tsv` whose `num_genes`,
  `num_cds`, `missing_genes` and `frameshift_flag` columns were empty for every sample. The
  five-sample invert regression panel ran green and reported nothing.

  The cause was two processes that were each both a QC step and a database writer, so the
  only way to reach the verdict was through the uploader. Both are split:

  - `annotation_stats.py` moves out of `PUSH_MTDNA_ANNOTATION_RESULTS` into a new
    `ANNOTATION_STATS`. The pusher keeps only `push_emma_annotation_results.py`.
  - `species_validation.py` no longer reads the nominal species from the DB (it comes from
    the samplesheet's `nominal_species_id`, already on meta) and no longer writes
    `lca_validation`. It emits a validation record that the new `PUSH_SPECIES_VALIDATION`
    upserts, following the existing `PUSH_ENA_VALIDATION_RESULTS` pattern.

  The DB-free chain -- `SPECIES_VALIDATION`, `ANNOTATION_STATS`, `EVALUATE_QC_CONDITIONS`,
  `QC_SUMMARY` -- becomes a new `MITOGENOME_QC` subworkflow that always runs, and
  `COMPILE_HELD_SAMPLES` moves out of the guard with it. `subworkflows/local/upload_results_mito`
  is now nothing but pushers. The rule the split enforces: a process either produces a verdict
  or writes to the database, never both.

  **Receipt gating.** `PUSH_LCA_BLAST_RESULTS` and `PUSH_MTDNA_ANNOTATION_RESULTS` used to
  inherit `gateOnAssemblyReceipt` transitively, by consuming channels derived from a
  `SPECIES_VALIDATION` whose input was gated. That inheritance is gone now that species
  validation is DB-free QC running upstream, so both take an explicit gate; without it the
  split would silently have reintroduced the foreign-key write-ordering race the helper
  exists to prevent. Every pusher is now gated at its own call site.

  **Ordering.** The receipt gate on the QC verdict itself is gone, because the verdict no
  longer writes anything. The one consumer that genuinely needs the ordering --
  `ENA_SUBMISSION_PREP`, whose `PREPARE_ENA_METADATA` reads `mean_depth` from
  `mitogenome_data` -- is now gated on `UPLOAD_RESULTS.out.assembly_receipts` at its call
  site in the main workflow.

  All seven database writers and all six tables they touch are unchanged.

  `nextflow_run_external.sh` passes `--skip_upload_results true` and will now get QC output
  it did not get before.

- The old `MITOGENOME_QC` subworkflow is renamed `ENA_SUBMISSION_PREP`, which is what it does:
  table2asn, flatfile conversion, candidate packaging and Webin validation. It runs downstream
  of the QC verdict and only ever sees samples the gate released, so it could not host the
  steps that produce that verdict. The name is now free for the subworkflow that does.


- Invertebrate runs no longer abort the entire session the moment the first `MITOS2` sample
  finishes annotating.

  `MITOS2.out.gff_proteins` gained a fourth element (`annotation/cds`, which the QC gate needs
  for its PCG ORF check) when the gate was hardened, but the filter feeding
  `ANNOTATION_QC_GATE` still destructured three. Groovy cannot spread a four-element list into
  a three-parameter closure, so it passed the whole list as one argument and threw a
  `MissingMethodException` from inside the `filter` operator. That is an operator failure, not
  a task failure, so `errorStrategy` did not apply: Nextflow aborted the session and killed
  every queued and running task. Both invertebrate test panels died that way with every task
  green (18/18 and 55/55, all `exit 0`), which is why the logs named no failing process.

  The closure now destructures all four elements and names the tuple shape at the call site.
  Every other 3+ element emit in the repo was audited against its consumers; this was the only
  mismatch.

  New `tests/coral_annotation_gate` runs the real `MITOGENOME_ANNOTATION` subworkflow fully
  stubbed over a coral, a non-coral invertebrate and a vertebrate, so every channel is bound
  and any future arity drift between a module's emit and the closure consuming it fails here
  instead of mid-run. The module tests could not catch this: they build the gate's input by
  hand and never exercise the wiring.

## v2.0.0 - [2026-09-03]

Second major release, and the first to carry the ENA submission path. Everything ENA-related is new since
v1.1.0, which shipped neither the params nor the tables: packaging, flatfile generation, Webin validation and
the whole `sql/` migration chain. The release also hardens post-assembly routing and the assembly summary,
makes the migration chain replayable end to end, and adds the annotation rescue and QC-gate work that closed
out the cycle.

Nothing here breaks a released version. The `Deprecated` section records parts of the ENA layer that were
built and then handed to a separate downstream pipeline within this same cycle, so no release ever carried the
params or tables it retires.

**Operator note.** The `sql/` migrations are applied to the live database as they land rather than at release
boundaries, so a database already carrying earlier migrations is affected regardless of when this tag lands.

The entries are grouped in two blocks. The first covers work that landed after the v2.0.0 changelog entry was
first drafted on 2026-08-18; the second is that original entry, kept intact. Both are part of this one release
-- the 2026-08-18 tag was never published.

### Work that landed after 2026-08-18

#### `Fixed`

- Gene and CDS extraction no longer labels every sample with the vertebrate genetic code.

  `bin/extract_genes_gff.py` and `bin/extract_cds_from_tbl.py` wrote a literal `[mgcode=2]`
  into every FASTA header they emitted, and neither module passed a code at all. They were
  the last step in the QC chain still ignoring `meta.genetic_code`. A code-4 coral therefore
  carried `[mgcode=4]` on its genome FASTA and `transl_table 4` in its `.tbl`, but `[mgcode=2]`
  on all 13 extracted genes -- and on the proteins, since `bin/translate_genes.py` copies the
  record description verbatim. Because the concatenated `*.genes.fa` is shipped in the ENA
  candidate package, single packages went out holding a code-4 genome record beside code-2
  gene records.

  Both scripts now take a required `--genetic-code`, validated against
  `orf_utils.SUPPORTED_CODES`, and both modules resolve it the way every other QC module does
  (`task.ext.code ?: meta.genetic_code`). There is deliberately no default: a wrong table
  reaching a submitted annotation is worse than a stopped task.

  Already-published assemblies keep the wrong tag until re-QCed -- see `docs/usage.md`.

- The run-level assembly summary no longer judges invertebrates against a vertebrate gene count.

  `--expected-gene-count` (37) is a vertebrate figure, applied run-wide. Cnidarians carry ~15
  genes because most of their tRNAs are nuclear-encoded, so finished coral mitogenomes were
  reported as `missing_genes` -- contradicting `annotation_stats.py`, which had already passed
  the same assembly on the protein-coding + rRNA core. `annotation_stats.py` now records the
  profile it used in a `completeness_profile` column and `mitogenome_assembly_summary.py`
  honours it. The 13-PCG check is true for these lineages too, so it still applies and still
  blocks.

- `annotation_stats.py` selects its completeness profile from the resolved genetic code rather
  than a hardcoded cnidarian class list. Code 2 is judged against the vertebrate 37-gene set
  and gene order; every other table against the conserved core. Keying on the class list meant
  a code-9 echinoderm or code-5 mollusc was judged against vertebrate gene order and failed for
  being what it is. `--class` remains as the fallback for callers with no code.

- `orf_utils.start_codons()` / `stop_codons()` raise on an unrecognised table instead of
  silently falling back. The old fallback returned `("TAA","TAG")` while its comment claimed
  the vertebrate code -- but code 2 also stops on `AGA`/`AGG`, so a mistyped table quietly got
  the wrong stop set. Stale argparse defaults that implied a table were removed from
  `translate_genes.py` (`--table`) and `mitos_to_emma.py` (`--code`).

#### `Added`

- `bin/audit_lca_db_coverage.py` and `bin/backfill_lca_uploads.py`: find published LCA results
  that never reached the database, and put them there without re-running the pipeline.

  The audit compares the row counts in each `mitogenomes/<OG>/<assembly>/lca/` against
  `blast_filtered_lca`, `lca` and `lca_raw_results`, reporting only assemblies whose published
  file has data and whose table has none -- a sample with no hits above threshold writes empty
  files and correctly has no rows. It exits non-zero when it finds a gap, so it can gate a run.

  The backfill rebuilds the `lca_combined` / `blast_combined` inputs `SPECIES_VALIDATION` would
  have produced, reusing that module's own concatenation helpers, and hands them to the same push
  scripts the pipeline uses. It reads from `mitogenomes/`, not `species_validation/`, so it still
  works on an archived run whose `work` and `species_validation` directories have been pruned. It
  pushes only the tables the audit found empty and never passes `--force`, so there is no
  superseded history to prune and re-running it once the gaps are filled finds nothing to do.

- An intron-split `cox1` is now rebuilt as `cox1_0`/`cox1_1` in `bin/coral_fix_bed.py`, the
  same reference-transfer treatment `nad5` already got.

  Most scleractinians carry the group I intron only in `nad5`, but some carry a second one in
  `cox1`, holding a LAGLIDADG homing endonuclease ORF (MITOS calls that ORF `lagli`, correctly).
  MITOS then annotates only one of cox1's two exons, and nothing downstream noticed: OG2361
  published an 873 nt `CO1` and a 291 aa `MT-CO1` (barely half a cox1), and the LCA was run on
  that fragment. It resolved to Scleractinia anyway, which is the dangerous part: half a CO1 still
  BLASTs to plausible neighbours, so the barcode looked fine.

  The reference's cox1 is normally single-exon, which is exactly what makes it a usable probe: it
  BLASTs onto an intron-split assembly as two subject blocks. `group_exons()` collapses HSPs that
  share a diagonal (so an exon running off the end of the contig and continuing past position 1
  stays one block) and treats a broken diagonal as the intron. The 5' end keeps the transferred
  start codon, the 3' end is walked in frame to the first stop, and the result must be a clean ORF
  before it replaces MITOS's row. On OG2361 this rebuilds cox1 as 17218-17928 + 730-1590, 1572 nt /
  523 aa, 99.0% identical to the *Favites abdita* NC_035879 COX1 over its full length, with the
  splice junction at the canonical anthozoan site (`...FWFFGH` | `PEVYIL...`). Across the 36
  batch-20 corals it is the only sample whose output changes.

- `bin/annotation_qc_gate.py` gained a `--min-co1-aa` (default 450) truncation check for cnidarians,
  mirroring the existing `--min-nd5-aa`. The truncated CO1 above did trip the gate's no-stop check,
  but only by luck of where the exon happened to end; the length now says so directly.

- Post-EMMA gene rescue for `ND4L` and `ATP8`: `modules/local/emma_gene_rescue_gate` +
  `modules/local/emma_gene_rescue`, driven by `bin/emma_rescue_gate.py` and
  `bin/rescue_emma_pcg.py`.

  EMMA's `rationalise_matches!` discards a short CDS when its computed circular overlap with a
  longer neighbour exceeds half the shorter feature's length, which routinely loses ND4L (against
  ND4) and ATP8 (against ATP6). The gene is in the assembly and EMMA even reports the match; only
  the annotation is short. Those bundles then fail `annotation_stats.py` (`passed=no`) and are held
  out of QC/ENA for a defect the assembly does not have.

  The gate reads the EMMA GFF and emits `FIX\t<targets>` only when the sole missing REF genes are
  ND4L and/or ATP8, both flanks of each are present, and every other REF gene is present and in
  order; anything else is `PASS\t-` and flows through untouched. The rescue then rebuilds each
  target from the flanking-gene coordinates EMMA already produced: define the intergenic window
  between the REF-order neighbours, tblastn a reference protein set into it to fix the reading
  frame, refine to a clean ORF (including EMMA's polyadenylation convention where the stop is
  completed by the poly-A tail), and write matching gene/mRNA/CDS lines into the `.gff`, the `.tbl`
  and the per-gene `cds/` and `proteins/` FASTAs. Every edit is guarded on BLAST identity and
  coverage, ORF cleanliness, and length against the matched reference; a target that fails any
  guard is left alone. Only the annotation bundle is swapped, so co1/12S/16S, BLAST and the LCA are
  untouched, and a rescue that recovers nothing re-emits EMMA's original bundle and the sample is
  held exactly as before.

  The reference set is `assets/panels/rescue_pcg_refs.faa` with `assets/panels/rescue_pcg_refs.manifest.tsv`,
  built by `bin/build_rescue_pcg_refs.py` from RefSeq mitochondrion CDS translations for
  Actinopterygii and Chondrichthyes plus a small tetrapod outgroup. The committed copy is the
  artifact the pipeline ships; the script is a stdlib-only refresher, not a runtime dependency.

- Post-EMMA tRNA rescue: `modules/local/trna_rescue_gate`, `modules/local/trna_scan` and
  `modules/local/trna_rescue`, driven by `bin/trna_rescue_gate.py` and `bin/rescue_trna.py`.

  EMMA's covariance model periodically misses a tRNA that is physically present on an otherwise
  complete, correctly ordered vertebrate mitogenome. The gate routes an assembly to the rescue only
  when every missing REF gene is a tRNA and the whole 13-PCG + 2-rRNA core is present and ordered,
  so each target's insertion gap is well defined. tRNAscan-SE 2.0 (vertebrate-mitochondrial model)
  is then run as a second, independent finder against the EMMA genome FASTA; the scan and the
  splicer are separate processes because tRNAscan's Perl container has no Python.

  A hit is spliced back only if it clears every guard: isotype *and* anticodon match the specific
  missing gene (which is what separates the two Leu and the two Ser isotypes), Infernal score above
  `--min-score`, length in range and intronless, midpoint inside the genomic gap between the
  target's nearest present neighbours (an origin-spanning gap is skipped), no more than
  `--max-overlap` bp of overlap with an existing feature, and exactly one surviving hit. Accepted
  hits are written as gene + tRNA lines into the `.gff` and `.tbl` the way EMMA writes its own.
  Failures are recorded per target in a status file and the script always exits 0.

- `annotation_trna_tolerance` (default 2): a vertebrate mitogenome that carries the whole conserved
  core (13 PCGs + both rRNAs) in the correct order but is short at most this many tRNAs now clears
  the QC/ENA gate instead of being held. That shortfall is an EMMA tRNA-model limitation rather than
  an assembly defect, and after the rescue above it is what is left over. `missing_genes` still
  lists every absent gene; the tolerated ones are additionally named in a new `trna_advisory` column
  (`sql/021_mitogenome_data_trna_advisory.sql`, plus the same field through
  `annotation_stats.py`, `evaluate_qc_conditions.py` and the push scripts) so a pass driven by the
  allowance stays queryable and auditable. `0` restores the old requirement of a complete 37-gene
  annotation.

- `bin/annotation_qc_gate.py` now runs a code-generic per-PCG ORF check for every invertebrate
  lineage: each of the 13 protein-coding genes is read from the spliced `annotation/cds/` FASTA
  MITOS wrote and checked for a valid start codon, a terminal stop, and internal stops against the
  sample's own translation table. This is exactly what table2asn enforces
  (`SEQ_FEAT.StartCodon`, `SEQ_FEAT.NoStop`, internal stop) and the check the gate previously
  lacked: a mis-placed boundary such as a MITOS ND1 off by three codons went straight through to
  table2asn and failed terminally there. The core-presence, ND5 and CO1 heuristics remain
  Anthozoa-specific and still run only for cnidarian (code 4) samples. The gate now takes
  `--cds` and a required `--genetic-code`.

- `bin/orf_utils.py`: one source of truth for the mitochondrial start/stop codon sets, covering
  every NCBI table the pipeline can route (2, 4, 5, 9, 13, 14, 21, 24, 33) with the per-table
  reasoning recorded. `bin/process_files.py` previously carried partial tables enumerating only
  code 2 and now reads them from here. Pure stdlib, so it imports in the gate's psycopg2 container.

- `bin/mito_gene_order.py`: the vertebrate `REF_GENES` order plus the tRNA / rRNA / PCG partitions
  and the tRNA anticodon and `/product` tables. `REF_GENES` had been copy-pasted into
  `annotation_stats.py` and the rescue scripts, each with a "keep this in sync" comment that had
  already drifted (three of them said "change both" or "change all three" while there were four
  copies). The gates and the QC step have to agree byte-for-byte on what "present and in order"
  means, so the list now lives in one place.

- A run-level `held_samples.tsv` (`modules/local/compile_held_samples`), one row per sample that
  did not reach submission-ready, with the cause. Two sources feed it: pre-QC holds, from a new
  machine-readable `held_reason` written by `evaluate_qc_conditions.py`
  (`species_not_in_blast`, `annotation_failed`, `assembly_anomaly:<type>`, `not_circular`), and
  table2asn quarantines with their blocking codes. Both sets were previously console-only
  `.view()` calls, so a run could report success with part of the batch quietly missing. The file
  is always emitted, header-only when nothing was held, so its presence is a reliable end-of-run
  signal rather than something that appears only on failure.

- `sql/020_ena_validation_attempts_og_num.sql` adds the generated `og_num` column to
  `ena_validation_attempts` as column 2, matching `sample`, `draft_genomes`, `lca`,
  `lca_raw_results`, `sequencing` and `mitogenome_data` (migration 018). PostgreSQL cannot insert a
  column at a position, so this is a table rebuild like 014 and 018. Migration 019 exists because an
  earlier rebuild silently reset a column for every row and nothing caught it until the corruption
  was found independently; this one verifies the copy row-for-row with a bidirectional `EXCEPT`
  diff over every carried column before the old table is dropped, so a mismatch raises inside the
  transaction and nothing is renamed.

- nf-test coverage for the new modules (`emma_gene_rescue`, `emma_gene_rescue_gate`, `trna_scan`,
  `trna_rescue`, `trna_rescue_gate`, `compile_held_samples`, `annotation_qc_gate`) and unit tests
  for `orf_utils`, `mitos_to_emma`, `coral_fix_bed`, `annotation_qc_gate`, `annotation_stats`,
  `evaluate_qc_conditions`, `species_validation`, and both rescue scripts and their gates.
- `qc_only_from_annotations.nf` now pushes the QC stage's own results to SQL through
  `UPLOAD_ENA_RESULTS` (`ena_validation_attempts` plus `lca_validation.validator_2`), instead of
  being entirely read-only. A sample QC'd through this entrypoint previously read as un-uploaded
  even though its QC had run.

  The scope is deliberately the QC stage alone: this entrypoint runs no assembly, no annotation
  and no LCA, so it writes no `mitogenome_data`, `blast_filtered_lca` or `lca` rows, and
  `prior_upload_status_files` is empty because the five pre-QC pushes never ran. Keeping
  `SPECIES_VALIDATION` off this path is the point of the narrow scope: under
  `--force_db_overwrite` that module overwrites `lca_validation.validated_species_name` and
  `validator`, which is exactly the hand-validated row this entrypoint exists to serve.
  `--skip_upload_results true` restores the old read-only behaviour. Pinned by
  `tests/qc_only_upload/main.nf.test`.


#### `Fixed`

- The QC-only entrypoint no longer gives every sample the vertebrate genetic code.

  `qc_only_from_annotations.nf` has no samplesheet, so it could not take `meta.genetic_code`
  the way the main pipeline does and assumed the run-level `--translation_table` for every
  sample instead -- code 2 unless the operator remembered otherwise. That is not confined to
  QC: `meta.genetic_code` becomes the `mgcode` in `GEN_FILES_TABLE2ASN` and the table in
  `TRANSLATE_GENES` and `FORMAT_FILES`, so re-QCing a coral assembly through this entrypoint
  rewrote an annotation made under code 4 as code 2 and submitted it that way.

  `VALIDATED_SPECIES_QUERY` now returns the taxonomic class alongside the validated species
  name -- same query, same round trip, no extra process -- and the entrypoint resolves the code
  from it. The class is looked up from the validated name first and only then from the sample's
  nominal one, since a hand-corrected species is exactly what this entrypoint exists for;
  matching is exact species then genus, deliberately not the fuzzy family/order tiers
  `bin/create_samplesheet.py` also has, which are loose enough to pick a reference but far too
  loose to pick a translation table.

  A class with no confirmed code falls back to `--translation_table` with a warning naming the
  sample and the class, rather than aborting as `prepare_samplesheet` does. The abort exists to
  stop a wrong table being baked into an annotation that is about to be made; these annotations
  already exist, so the useful thing is to run and say which samples need the flag.

- The class -> genetic code lookup moves out of `prepare_samplesheet` into
  `lib/MitoGeneticCode.groovy`, so the main pipeline and the QC-only entrypoint resolve one
  table instead of two. Only the lookup is shared: what to do with an unmapped class stays with
  each caller, because they genuinely differ. Pinned by `tests/mito_genetic_code`.

- A row PostgreSQL rejects no longer costs a sample its entire LCA upload.

  `bin/push_lca_blast_results.py` and `bin/push_lca_raw_results.py` wrapped each row's INSERT in
  its own `try`/`except`, which looks like per-row isolation but is not: PostgreSQL aborts the
  whole transaction on any error, so every later row failed with `current transaction is aborted`
  and the commit was downgraded to a rollback. The scripts then printed a tally of "succeeded"
  rows and exited 0, so the loss left no trace -- the Nextflow task showed COMPLETED and its
  published `.upload.txt` ended in a tick, while the table had nothing.

  Each row now runs inside a savepoint (`bin/pg_row_guard.py`), so a rejected row costs that row
  alone and the rest of the batch commits. Both scripts exit non-zero when any row failed, and
  the tick is reserved for an upload that landed in full; a partial upload says
  `finished with errors`. `tests/unit/test_push_lca_row_isolation.py` reproduces PostgreSQL's
  abort semantics against a fake cursor, so the pre-fix behaviour fails it.

  This cost 87 assemblies their LCA rows across `batch-12` .. `batch-20`, in three classes: a
  `;`-joined `staxids` against `blast_filtered_lca.taxon_id` (integer), a confidence value of
  ~1e-163 against `lca.top_confidence_score` / `lca_raw_results.confidence_score` (`real`, floor
  ~1.18e-38), and an insert that arrived before its `mitogenome_data` parent and tripped the
  foreign key. Migration 022 widens the three columns; `bin/backfill_lca_uploads.py` restored the
  rows from the published output, with no pipeline re-run.

- A MitoHiFi run that crashed after writing its assembly no longer disappears from the pipeline.

  `MITOHIFI_CHECK_CIRCULARITY` inner-joins MitoHiFi's `coverage_mapping` output, which is
  `optional: true`. A run that died between the final FASTA and the coverage step therefore
  produced no circularity evidence, and the `ch_circ_verdict` join then discarded the sample
  outright: no annotation, no QC summary, no `held_samples.tsv` row, no database row. The run
  still reported success, because the module wrapper runs under `set +e` so the task exits 0
  whatever MitoHiFi did.

  OG2133 (*Benthalbella* sp.) is the case in point. It assembled a circular 21,236 bp contig in
  batch 18 and then vanished: of 69 samples it was the only one with no `mtdna.upload.txt`. The
  crash was `KeyError: 'product'` in MitoHiFi's own `getGenesList.py`, parsing the reference it
  had been given -- AP012968.1 carries a tRNA annotated `/note="tRNA-undetermined"` with no
  `/product`, which MitoHiFi dereferences unconditionally. Any NCBI record with a
  qualifier-less feature reproduces it.

  `ch_mitohifi_fasta_branched` gains a `partial` arm for a non-empty FASTA whose run left a
  traceback in the command log, and that arm is now a complement in every channel that
  establishes totality (`ch_circ_verdict`, `ch_circularity_evidence`, `ch_assembly_log`,
  the summary inputs), so the sample reaches annotation and is held visibly at QC instead of
  being dropped. It keeps its published reference-relevance diagnostic but stays out of the
  Oatk fallback routing, which needs circularity evidence it never produced.

  A `partial` sample carries MitoHiFi's own `was_circular` rather than a blanket unknown.
  `bin/check_circularity.py` computes `verdict = (mitohifi_circ is True) || read_span || hifiasm`,
  a monotone OR, so the check can only ever flip `False -> True` -- across the 27 checked
  assemblies in batch 18, `True -> True` 19 times and `False -> True` 8 times, never the
  reverse. A `True` is therefore exactly the verdict the check would have reached and is
  trusted; a `False` degrades to unknown and never to `false`, because that is the
  terminal-overlap false negative the check module exists to repair and it is unresolvable
  without the coverage mapping.

- Every test in `tests/assembly_routing` was failing before it ran.

  `conf/base.config` reads `params.mitogenome_depth_max_forks` at parse time, added with the
  coverage fork-bounding, but the suite's `nextflow.config` pre-declared only `params.outdir`.
  All five oatk-fallback tests died with `Unknown config attribute` before reaching a workflow.
  The param is now pre-declared alongside `outdir`, with a note to keep the list in step.

- `MITOGENOME_COVERAGE`, `OATK`, `LCA` and `SPECIES_VALIDATION` now retry a walltime kill instead
  of silently dropping the sample.

  All four carried a narrowed `errorStrategy` treating only `exitStatus == 137` as transient, so a
  Slurm timeout (140/143) fell straight through to `'ignore'`. This is the same defect fixed for
  `GETORGANELLE_RESEED`, where it cost OG28.ilmn.231024 its reseed assembly in batch-02; the
  narrow overrides are removed rather than widened, so all four inherit the global default
  (`maxRetries = 2`, transient `(130..145) + 104 + 247`, `'ignore'` on exhaustion) and there is one
  copy of the policy to maintain. Both original intents survive: `137` is still covered because it
  sits inside `130..145`, and depth is still never a reason to fail a sample because the inherited
  default also ends in `'ignore'`, not `'terminate'`.

  This was live risk, not theory. Batch-19 had three coverage tasks at 51-64% with more than five
  hours of work left against a 3h20m remaining walltime; under the old strategy each would have
  been ignored on timeout and published an empty depth placeholder.

- `MITOGENOME_COVERAGE` concurrency is now bounded by `params.mitogenome_depth_max_forks`
  (default 4).

  Each task streams an entire WGS library (~150 GB gzipped, ~2 billion reads) past a 33 kb doubled
  index through `minimap2 -t N | awk | python`, so a batch is limited by filesystem read bandwidth
  and by that single-threaded awk, not by cores. Nothing bounded submission, so batch-19 launched
  13 within four minutes and every one of them degraded roughly 15x, from 476k to 22k reads/s, with
  minimap2's CPU multiplier falling from 2.4x to 0.70x. The same samples, on the same nodes, with a
  byte-identical command, had held 290-430k reads/s to completion at 6 concurrent: OG2288 took
  1h46m then and had not finished after 4h40m under contention.

- `bin/build_source_modifiers.py` no longer invents a hemisphere for a latitude the sample
  table records without one.

  The sample table stores latitudes as unsigned magnitudes: a Ningaloo sample at 22.03 S is
  `latitude_collection = 22.03`. `parse_coordinate()` defaulted an unsigned value to the
  hemisphere implied by its sign, so every one of them came out N. The result is well formed,
  so `valid_lat_lon()` could not catch it; table2asn accepted the shape and then rejected the
  value as `SEQ_DESCR.LatLonValue` ("Latitude should be set to S") because the coordinate
  contradicts the country. In batch-19 that quarantined 15 Western Australian assemblies,
  OG2279 through OG2294, the entire 260114 ilmn plate. The two parse paths also disagreed
  about the default: `fallback_parse_latlon()` assumed S and E, which would have been right
  here, but `smart_split_latlon()` splits a bare `"22.03 113.891"` pair successfully so the
  fallback never ran.

  The hemisphere is now only ever read, never inferred from a sign that is not there.
  `parse_coordinate()` takes the axis (so a longitude written `23.43 S` is rejected rather
  than emitted as a second latitude) and returns `(magnitude, hemisphere-or-None)`, where
  None means the value states no hemisphere, which is a different thing from an unparsable
  one. A new `format_lat_lon()` resolves that None against `COUNTRY_HEMISPHERE` in
  `bin/geo_loc_name_utils.py`, keyed on the *resolved* geo_loc_name so the alias table stays
  the single place country spellings are corrected. Where the country cannot settle it, the
  modifier is omitted and the reason named: the country straddles that axis (Indonesia,
  Brazil, Ghana, Kiribati), or is not a mapped country at all. A stated hemisphere the
  country contradicts is likewise dropped rather than corrected, since which of the two is
  wrong cannot be told from here. That turns a hard `FAIL_TABLE2ASN` into an omitted
  optional qualifier, so one bad row can no longer quarantine an otherwise clean assembly.

  Country resolution consequently moved above the coordinate block in `main()`, and the old
  post-hoc `unknown`/`valid_lat_lon` blanking passes are subsumed by `format_lat_lon`.

  This makes the pipeline resilient to an unsigned latitude; it does not make the stored
  value right. Preserving the sign in the spreadsheet ingest remains the durable fix.

- `bin/rotate_to_cox1.py` no longer re-origins inside cox1.

  The script's premise was that cox1 is "a conserved, single-exon gene corals always carry" and "is
  not the wrapping feature". For a coral with an intron-split cox1 both halves of that are false. It
  anchored on the best-scoring tblastn HSP, which is usually the 3' exon, then back-extrapolated
  `3*(qstart-1)` to reach the N-terminus, landing in the middle of the intron. On OG2361 that put
  position 1 at offset 568, cutting cox1 across the origin: precisely the failure the rotation exists
  to prevent, moved from nad5 to cox1.

  The anchor is now the lowest-`qstart` HSP of the best-scoring query (cox1's true 5' exon), and a
  new `MAX_EXTRAPOLATE_AA` refuses to rotate at all when even that HSP starts too far into the
  protein to locate the 5' end. OG2361 now re-origins at 17768, cox1's actual start, leaving the
  gene linear at 19-729 + 1808-2668. It is the only batch-20 sample whose rotation changes.

- The run now aborts when a sample's taxonomic `class` is unresolved, and `class`/`family`/`order`
  resolve from the NCBI taxdump when the `species` table has no match for the sample's nominal
  name.

  `unknown` was never a neutral value. `is_invertebrate()` maps it to `false`, so an unresolved
  sample silently took the vertebrate path through every chooser downstream: genetic code 2
  instead of 4/9, EMMA instead of MITOS2 (and no `ROTATE_ORIGIN`), the curated fish BLAST DB
  instead of `nt`, the 82.0 rather than 88.0 reference-relevance threshold, and the vertebrate
  22-tRNA completeness expectation pushed to the database as a QC verdict. A coral run that way
  produces results that look normal and are wrong. `PREPARE_SAMPLESHEET` now checks the whole
  parsed sheet before any work is dispatched and fails with every offending sample named;
  `--allow_unknown_taxonomy` downgrades it to a warning for a deliberate one-off. A blank
  `family`/`order` with a known class only warns, since it degrades the reference-divergence
  tiering rather than flipping invert/vertebrate routing.

  The samplesheet is still written and published before the abort — it is the artefact to
  correct — and now lands in `<outdir>/samplesheet/` (the existing `withName: CREATE_SAMPLESHEET`
  selector is a full match and never covered `CREATE_SAMPLESHEET_ENRICHED`, which had been
  publishing to `<outdir>/create/`).

- `bin/taxdump_lineage.py`, a stdlib-only `nodes.dmp`/`names.dmp` resolver, wired into
  `CREATE_SAMPLESHEET_ENRICHED` as a fourth input from the already-cached
  `DOWNLOAD_TAXONKIT_DB` (aliased, `storeDir` makes it a cache hit). The species table only
  holds curated taxa, so it misses valid names outright: `OG85` / *Epinephelides armatus*
  resolved to `unknown` from the database and now resolves to Actinopteri / Serranidae /
  Perciformes. Lookup is by scientific name, tries the binomial then the genus, refuses
  cross-kingdom homonyms, and skips NCBI's open-nomenclature placeholder nodes
  (`Acropora sp.`) in favour of the genus. The curated database still wins wherever it has an
  answer; the taxdump only fills blanks. Per-sample provenance is written to
  `taxonomy_resolution.tsv` next to the samplesheet rather than to a new samplesheet column,
  which would have changed every `meta` map's shape and silently emptied joins on resume.

- `PUSH_QC_VALIDATOR` (`bin/push_qc_validator.py`, `modules/local/upload_results/qc_validator/`),
  which records the pipeline as the *second* species-ID validator in `lca_validation`. A
  mitogenome needs two validators signed off before it is OK to submit; the pipeline already
  filled the first (`validator = 'nf-core'`, from the LCA/BLAST check in `SPECIES_VALIDATION`)
  and the second was a manual step that bottlenecked every batch. Any sample reaching
  `submission_ready = true` — table2asn PASS, EMBL flat-file conversion PASS, and the webin-cli
  format check PASS — now gets `validator_2 = 'QCd-nf-core'` written automatically from
  `UPLOAD_ENA_RESULTS`, off the same `<full_seqid>.ena_validation_result.tsv` that
  `PUSH_ENA_VALIDATION_RESULTS` consumes.

  The write is deliberately one-way. The `UPDATE` is guarded on `validator_2` being NULL or
  blank and there is **no `--force`**: the column exists to record that a second, independent
  validator looked at the sample, so clobbering a human's sign-off is never the right move. It
  also never INSERTs — with no `lca_validation` row there is no first validator either, and a
  lone `validator_2` would mean nothing — so a missing row is logged and skipped. Failures are
  non-fatal (`set +e` / `exit 0` / `UPLOAD_EXIT=`), like every other push module. Note that
  `submission_ready` is always `false` under `--ena_webin_validate false`, so the second
  validator is not recorded in runs with the format check switched off.

- `assets/placeholders/ena_not_run/`, the six NOT_RUN placeholder files the fixed-totality ENA record grouping
  resolves with `checkIfExists: true`. They were referenced but never created, which aborted every
  run during workflow construction. Their contents are inert: `collate_ena_validation.py` reads
  only the real status-file suffixes and infers NOT_RUN from a stage's absence, so the placeholders
  exist purely to make each sample's contributor count fixed and known.

- `precomputed_mitogenome_assembly_fasta_oatk` and `precomputed_mitogenome_assembly_log_oatk` are
  declared in `nextflow_schema.json`. They have worked from `nextflow.config` all along, but being
  undeclared meant the schema validator flagged them as invalid input on every run.

- `ena_submissions`, a submission ledger, plus the `ena_submission_status` view
  (`sql/015_ena_submissions.sql`). Until now the database could say a flatfile cleared every gate
  but not what happened to it afterwards: whether it was submitted, what ENA returned, or which
  BioSample the submission was registered against. That lived only as receipt XML under `receipts/<OG>/` in the
  downstream submitter plus `mitogenome_data.genbank_accession`, which is keyed without the
  annotation version and named for an archive that does not mint ERZ accessions. The ledger is
  deliberately a separate table rather than columns on `ena_validation_attempts`: that table is
  upserted with `DO UPDATE SET` across every non-key column, so an accession stored there would be
  erased by the next validation rerun. Keyed on `(full_seqid, webin_mode)` so a re-annotation is a
  new row and a test-service dry run cannot overwrite the production record; identity columns are
  `GENERATED ALWAYS` from `full_seqid` so they cannot drift. It holds submission status
  (`NOT_SUBMITTED`/`SUBMITTED`/`ACCESSION_ASSIGNED`/`FAILED`, with constraints requiring the
  evidence each status claims), timing, receipt path and digest, the ENA accessions (`ERZ`, `GCA`,
  sequence, `ERS`, `PRJEB`), the BioSample actually used and its source, the locus tag prefix, and
  the run accessions. **This pipeline never writes or reads it** — validation must not depend on
  submission state — so the writer is the downstream submitter; `docs/ena_submission_handoff.md`
  carries the contract and the idempotent-upsert pattern.

- The mitochondrial translation table is now resolved once and asserted, instead of being defaulted
  to 2 or 4 at each point of use.

  `mitoGeneticCode()` in `subworkflows/local/prepare_samplesheet` gained the sample id (so its
  errors name the sample), accepts an explicit per-sample `genetic_code` samplesheet column that
  wins over the class map, and now raises for an `invertebrates=true` sample whose class it does not
  know rather than silently returning the vertebrate default. As more invertebrate lineages are
  added (bivalves are code 5, for instance) a wrong default would mistranslate an entire annotation
  and fail table2asn terminally. `MITOGENOME_ANNOTATION` then asserts, before the EMMA/MITOS2
  branch, that every sample carries an integer `meta.genetic_code` in the supported set, mirroring
  the existing `--mitos_refdb` / `--nt_blast_db` asserts.

  With the value guaranteed upstream, the `meta.genetic_code ?: 2` and `?: 4` fallbacks scattered
  through `mitos2`, `annotation_qc_gate`, `coral_annotation_fix`, `translate_genes`,
  `gen_files_table2asn` and `format_files` are gone; each now reads `task.ext.code ?:
  meta.genetic_code`, so a missing code fails loudly at the assert instead of quietly annotating
  under one table and validating under another.

- `bin/mitos_to_emma.py` now carries the tRNA anticodon through to the `/product` string, so a
  MITOS2 `trnW(tca)` becomes `tRNA-Trp(UCA)` exactly as EMMA writes it. `map_gene_name()` returns
  the anticodon alongside the name, type and fragment label, and the product is built from the
  shared table in `bin/mito_gene_order.py` so the two annotators cannot emit different strings for
  the same feature. An unrecognised suffix or anticodon falls back to a bare gene-name product; a
  cosmetic field should never fail the run.
- A walltime kill on `GETORGANELLE_RESEED` is retried again. Its `withName` block in
  `conf/base.config` overrode the global error strategy and narrowed "transient" to exit 137
  (OOM) alone, so exit 140 -- what Slurm returns when a task exceeds its time allocation --
  fell straight through to `ignore` on attempt 1, with `maxRetries = 2` never spent.
  `GETORGANELLE_RESEED` was the only process carrying that override; `GETORGANELLE_FROMREADS`
  next to it has always inherited the global policy, which counts 130..145 as transient.

  It cost a real assembly: in `batch-02`, `OG28.ilmn.231024` timed out at 16 h and was ignored,
  leaving `mitogenomes/OG28/OG28.ilmn.231024.getorg1770reseed/mtdna/` holding nothing but
  `reference_seed/` -- no FASTA, no GFA, and therefore no downstream annotation for that library.
  The override is deleted rather than widened, so the two GetOrganelle processes now resolve to
  one retry policy and 137 stays covered because it sits inside 130..145. Note the second attempt
  gains only 8 hours (`Math.min(24, 16 * task.attempt).h`, capped by `max_time = 24.h`): samples
  that exceed 16 h are usually not converging rather than running slightly long, so expect the
  retry to buy a verdict rather than an assembly.

- `WEBIN_VALIDATE` bounds each `webin-cli` call with `timeout` and retries transient failures
  in-script, controlled by `--webin_validate_timeout_seconds` (default 900) and
  `--webin_validate_max_attempts` (default 3). The call was previously unbounded, so a hung
  ENA request could only end when Slurm killed the whole task: `OG16` in `batch-02` spent its
  entire 4 h allocation inside a call that normally returns in about 5 seconds, and the
  automatic task retry then passed in 4.

  Only timeouts and infrastructure failures are retried, with a 30 s/60 s backoff; a
  `FAIL_WEBIN` verdict is deterministic and breaks out immediately rather than re-running a
  rejected flatfile. `webin_output` is cleared between attempts, since the classifier greps it
  and a report left by an earlier attempt would misclassify a later one. The process still
  ends `exit 0` and the `.webin_status.tsv` schema is unchanged -- a failing task would drop
  the non-optional `manifest`/`status`/`log`/`reports` emits and erase the sample from the ENA
  collation instead of recording it as failed -- so the only new value is the `webin_timeout`
  reason, and the attempt count goes to the log.

  The biocontainer ships **busybox** `timeout`, not GNU coreutils: it takes positional seconds
  with no `--signal`/`--kill-after`, and it exits **143** on timeout rather than GNU's 124. The
  module matches both codes so the `conda` path, which does supply GNU `timeout`, behaves the
  same.
- `BUILD_SOURCE_MODIFIERS`'s stub emitted files under names that did not match its `output:`
  block (`output/bankit_metadata.csv`, `src_files/dummy.src`). The optional `src_file` emit
  therefore produced nothing, the join into `GEN_FILES_TABLE2ASN` came out empty, and every
  `-stub` run of `MITOGENOME_QC` was silently truncated before `table2asn` -- the workflow still
  reported success. The stub now writes `${meta.id}.bankit_metadata*.csv` and
  `${meta.mt_assembly_prefix}.stub.src`, matching the declared outputs.

#### `Changed`

- `MITOGENOME_COVERAGE` is sized from measurement: `cpus` 12 -> 4 and `memory` 16 GB -> 4 GB per
  attempt. The Aug-25 batch-19 trace for OG2288 records `%cpu=2228` (2.2 cores) and
  `peak_rss=484604` (473 MB), and threads past the single-threaded awk consumer cannot be used.
  Memory is the directive that actually shrinks the allocation, since Slurm drives the granted core
  count up to satisfy the memory request.

  Neither change affects `-resume`. Nextflow's task hash covers session id, process name, script
  source, container fingerprint, conda env, module, arch, bin dirs, inputs and attempt number;
  resource directives and retry policy are not hashed. The `memory_hints.json` preamble in
  `conf/base.config` said otherwise and has been corrected: the reason a resumed retry misses its
  cached success is the attempt salt, not the resolved memory. The workaround it documents is still
  necessary and unchanged.

- `apply_ena_migrations.py` now applies each migration exactly once, recording it in a
  `public.schema_migrations` ledger (filename, sha256, applied_at) and skipping anything already
  there. It used to re-execute the whole list on every invocation. Every migration is written to be
  idempotent so that was survivable, but not harmless: `003`'s label-only backfill has no date guard,
  so each run stamped `legacy_getorg_kmer` onto whatever rows happened to have a NULL `depth_method`
  at the time. That hit the same five current-era rows twice and had to be undone by hand both times.
  Any future migration inherited the same blast radius — a one-column change quietly mutating
  unrelated data.

  `--baseline` records the listed migrations as applied without running them, to onboard a database
  they were already applied to, and refuses once the ledger is non-empty. A migration whose file has
  changed since it was applied is now an error rather than a silent re-run, since the database no
  longer matches the file; `--force` downgrades it to a warning. `--check-only` reports the pending
  list.

- `sql/018_mitogenome_data_og_num_first.sql` puts `og_num` back as column 1 of `mitogenome_data`.
  `016` had to drop and re-add the column to make it generated, which moved it to position 91, and
  the position was deliberate. PostgreSQL cannot reorder columns, so `018` rebuilds the table with
  the columns declared in the wanted order: it copies the rows, refuses to commit a short copy,
  recreates the index, re-points the three inbound foreign keys from `lca`, `lca_raw_results` and
  `lca_validation`, re-grants `SELECT` to `readonly` and restores the column comments. The `_1`
  suffixes on the constraint and index names are preserved deliberately — the unsuffixed names still
  belong to the `mitogenome_data_SS260818` snapshot. The audit now also requires `og_num` to be
  column 1, since a rebuild that moves it is the same accident that lost the generation expression
  the first time.

- `sql/017_lca_content_addressed_rows.sql` writes down a schema change that had only ever been
  applied by hand: the `content_hash` and `taxon_rank_db` columns on `lca` and `lca_raw_results`,
  the `lca_set_content_hash()` trigger function and its four triggers, both unique keys swapped
  from `lca_run_date` to `content_hash`, and both foreign keys renamed off the ambiguous shared
  `fk_mitogenome`. The code half shipped in `55321df`, the same commit that added migrations 014
  and 015, so the ENA side was captured and this side was not.

  It mattered because `push_lca_blast_results.py` and `push_lca_raw_results.py` name
  `lca_content_unique` and `lca_raw_results_content_unique` as `ON CONFLICT` targets, so a database
  rebuilt from `sql/` failed on the first LCA upload of a run. Every step is guarded, so applying
  it to the existing database is a verified no-op: no column, constraint, index, trigger or
  `content_hash` value changed. `apply_ena_migrations.py` now audits for it.

- `sql/016_mitogenome_data_og_num_generated.sql` restores `mitogenome_data.og_num` as a stored
  generated column, `(SUBSTRING(og_id FROM 3))::integer`, matching `sample`, `draft_genomes`,
  `lca`, `lca_raw_results` and the `mitogenome_data_SS260818` snapshot, all of which are 100%
  populated. The live table had become a plain nullable column with no default, no generation
  expression and no trigger, and nothing writes it: neither `push_mtdna_assm_results.py` nor
  `push_emma_annotation_results.py` lists it in their upserts, so it was populated on 28 of 289
  rows, all of them pre-dating the loss. `mitogenome_submission_view` selects `og_num` and was
  reading NULL for most of the pipeline's own output as a result.

  PostgreSQL 14 cannot convert a column in place, so the migration drops and re-adds it, guarded
  so a re-run is a no-op. Safe here because `og_num` is in no key, index or foreign key, and no
  view depends on the live table's copy. Side effect: `og_num` moves to the last column position,
  where it already sits in `sample` and `draft_genomes`. `apply_ena_migrations.py` now audits the
  generation expression, so a table that loses it again fails the post-migration check instead of
  quietly filling with NULLs.

- The SQL upload summary gained a `Validator 2` column (`qc_validator` in
  `upload_results_summary.tsv`), reporting per sample whether `validator_2` was set, was
  preserved because someone already signed off, was skipped as not submission-ready, or found no
  `lca_validation` row. `compile_upload_report.py`'s `ENA_STEP` constant became the
  `SEQID_KEYED_STEPS` set, since two steps are now named after `meta.full_seqid` and both need
  the annotation token reconciled away before grouping — without that the new file would have
  produced a second, half-empty row per assembly.

- `<full_seqid>.manifest.txt` is gone and `<full_seqid>.package_metadata.json` is the whole ENA
  handoff. The standalone manifest could never be correct: `STUDY` is the BioProject the
  submission pipeline registers, `SAMPLE` is the BioSample it registers, and `RUN_REF` belongs to
  the raw-read submissions, so three of its twelve keys were values this pipeline does not have.
  Every published package was in fact the four-key `# BLOCKED` stub, and webin-cli had never
  checked one: genome-context validation resolves `SAMPLE` before anything else, so the only
  webin-checked artefact is, and always was, the sequence-context `<full_seqid>.webin_manifest.txt`
  under `ena/validation/flatfile/`. Package metadata goes to `schema_version` 3: a `manifest` block
  renders the nine keys this pipeline can fill, under the names Webin uses, so the submitter reads
  one file and adds its own three; a `specimen` block carries the flatfile's source-feature facts
  (organism, isolate, tissue, geo, date, lat_lon where present), read back off the packaged EMBL so
  the two can never disagree; `study` becomes `validation_study`, which is what it always was, the
  study sequence-context validation ran against and never a submission target; and
  `biosample_accession`, `biosample_source` and `run_accessions` are dropped along with
  `ena_package.py`'s `--biosample` and `--run-accession`. A field that fails validation is now
  omitted from the `manifest` block instead of collapsing the whole manifest to a stub.
  `prepare_ena_metadata.py` goes to `schema_version` 2 and no longer queries
  `sample.ncbi_biosample_id`. `package_digest` and `checksums.sha256` change for every package,
  which is correct: the contents changed. `refresh_package` carries an older package forward rather
  than leaving a hybrid, and deletes a legacy manifest file if it finds one.
  `docs/ena_submission_handoff.md` is rewritten against this shape. It was written in the same
  commit that changed the behaviour underneath it, so it still described the retired
  per-technology child studies, a manifest to read rather than build, and a BioSample field to
  check; its `in the flatfile?` column also read as permission to omit manifest keys, which it
  never was, since webin-cli checks the manifest and the flatfile independently. Every example
  in it is now copied from a real package rather than composed.

- The ENA study is one run-wide `--ena_study`, replacing `--ena_study_hifi`/`--ena_study_hic`/
  `--ena_study_ilmn`. The per-technology child studies existed because the umbrella PRJEB110568
  cannot receive data, so each technology needed its own study resolved from the candidate's
  technology. The study is now created by the upstream pipeline and handed to this one as a single
  accession, so technology no longer selects between studies.
  `subworkflows/local/utils_ena_targets` collapses to `validateEnaStudy`/`enaStudyAnnotate`; there
  is still no default, so a missing or non-`PRJEB` accession stops the run rather than validating
  candidates against whatever the run fell back to. The strict `sequencing_type` check that lived
  there is not lost: `platform_for_tech` in `bin/prepare_ena_metadata.py` already raises on any
  unsupported technology, so the failure simply moves from channel-build time to task runtime.

- ENA validation records are keyed on `full_seqid`, not the assembly prefix. What the pipeline validates is a
  flatfile built from **one annotation of one assembly** (`OG82.ilmn.240313.getorg1770.emma102`), and the EMBL,
  chromosome list, manifest, package and metadata JSON have always been named that way; only the validation record
  used the 4-field `mt_assembly_prefix`. So a re-annotation overwrote the record for the annotation validated
  before it (the upsert conflict target was `(assembly_prefix, ena_study, validation_attempt)`), both annotations
  of one assembly collided in the `groupTuple` key that releases records, and no row could say which annotation
  earned `submission_ready`. `assembly_prefix` is replaced by `full_seqid` in the record and in
  `ena_validation_attempts`, the identity split gains a fifth field, and the annotation version sits in its own
  `annotation` column after `code`. `<full_seqid>.ena_validation_result.tsv`, the Webin manifest/status/log and the
  table2asn status/findings now share the stem the flatfile already had; publish directories are unchanged.
  `ena.nf` takes an optional `full_seqid` CSV column, defaulting to `mt_assembly_prefix` (annotation then empty).
- `sql/014_ena_validation_attempts_full_seqid.sql` rebuilds `ena_validation_attempts` for the new column order,
  moves `ena_validation_attempts_key_idx` and `ena_validation_latest` onto `full_seqid`, and adds `annotation` to
  the identity index. Rows written before it are carried over with `full_seqid = assembly_prefix` and a NULL
  `annotation` (their annotation is not recoverable from the table) and are superseded the next time those
  assemblies are validated. The `assembly_prefix`-dependent statements in `001`, `002`, `011` and `012` are now
  guarded on that column still existing, so `bin/apply_ena_migrations.py` can keep replaying the whole chain.

### Work from the first half of the cycle (internally tagged 2026-08-18)

This is the v2.0.0 entry as first drafted, kept intact.

**Operator note.** The `sql/` migrations are applied to the live database as they land rather than at release
boundaries, so a database already carrying `001`-`009` is affected by `010`-`012` regardless. `010` archives
`ena_locus_registry` and `ena_candidate_loci` before dropping them; `012` drops `ena_candidate_packages`,
`ena_submission_selections` and the `ena_submission_queue` view with no archive, so dump those first if their
contents matter.

#### `Added`

- Uniform, cross-platform mitogenome read depth: `MITOGENOME_COVERAGE` + `bin/mito_depth.py`. `mitogenome_data.avg_coverage`
  previously held three different quantities depending on which assembler produced the row, so comparing it across
  platforms was always wrong: GetOrganelle wrote **k-mer** coverage (~0.2x true depth, and measured over the reduced read
  set its `--reduce-reads-for-coverage` default selects), MitoHiFi wrote per-base depth of *only* the reads recruited by
  mapping to a related-species reference (so a divergent reference silently depressed it, the same failure mode the
  divergence guard targets), and Oatk wrote nothing at all. One definition now replaces all three: mean per-base depth of
  the sample's own reads remapped to the assembly that actually reaches annotation. Remapping to **self** rather than to a
  reference is the point, since a genuinely divergent mitogenome is then measured just as accurately as a well-referenced
  one. Circular molecules are mapped to a head-to-tail doubled reference and folded (`d[p] += d[p+L]`), removing the false
  depth dip at both ends of the linearised molecule: on a synthetic 16.5 kb circle at a known 145.45x this recovers exactly
  145.45x folded versus 129.2x unfolded, and drops the CV from 0.22 to 0.08. NUMTs are rejected by **gap-compressed**
  identity rather than raw `NM`/aligned-length, so a read spanning a real control-region indel survives (`74M30D46M`,
  `NM:i:30`: raw 0.75, gap-compressed 0.992) instead of punching a depth hole in the D-loop. MAPQ is deliberately ignored,
  since a doubled reference gives every read two equally good placements. Runs once per sample on the `SANITISE_FASTA`
  output, so nothing is spent remapping assemblies that failed, fell below the length floor, or lost to a reseed. Results
  land in `<prefix>.mito_depth.tsv`, the assembly summary, MultiQC, and the new `mean_depth` / `depth_cv` / `breadth_*` /
  `mito_read_fraction` / `depth_method` columns in SQL (`sql/003_mitogenome_data_uniform_depth.sql`). Controlled by
  `--skip_mitogenome_depth`, `--mitogenome_depth_min_identity_{sr,hifi}`,
  `--mitogenome_depth_min_aligned_frac_{sr,hifi}` and `--mitogenome_depth_subsample_fraction`.
- `high_coverage_variability` is now advisory on a complete-core assembly, joining `low_mean_coverage` in
  `ADVISORY_WHEN_COMPLETE`. It was previously blocking only because it could barely fire: `coverage_cv` was populated for
  MitoHiFi alone, and even there it was measured against a *linear* reference, so 15-20 kb HiFi reads produced a
  triangular depth profile whose CV was mostly an artefact of linearisation. Now that the folded remap measures it
  properly for all three assemblers, leaving it blocking would mean measuring coverage *better* caused finished
  mitogenomes that previously passed to start failing. It stays visible in `manual_review_reason`, and still blocks on
  anything that is not a complete-core assembly. `parse_coverage` also now scans sources in an explicit order
  (`.mito_depth.tsv`, then `.coverage.tsv`, then `.contigs_stats.with_coverage.tsv` read from its `final_mitogenome` row);
  the previous predicate matched `.contigs_stats.with_coverage.tsv` through a loose substring test, so which file won was
  non-deterministic.
- Reference re-selection: `REFERENCE_CANDIDATES` + `REFERENCE_RANK` (`bin/rank_reference_candidates.py`).
  `findMitoReference` stops at the *first* complete mitogenome it meets walking up the sample's NCBI lineage, so a taxon
  with no congeneric record gets an arbitrary member of whatever rank the walk reached, and nothing in the pipeline ever
  looked for a better one. When the pre-assembly divergence guard reports a non-congeneric reference, these modules now
  fetch `--n_reference_candidates` candidates (the same `findMitoReference.py`, asked for more of them with `-n`) and keep
  the one a subsample of the sample's *own reads* maps to best — the quantity that actually determines whether MitoHiFi's
  reference-guided read recruitment succeeds, and one that is measurable before assembly. Scored by aligned read bases per
  reference base so a longer reference cannot win on length; ties keep the taxonomically closest candidate, so the choice
  can only improve on the previous behaviour. Wired into both the MitoHiFi route and the GetOrganelle vertebrate reseed
  (where the reference drives both the seed and the custom gene DB). A published `<prefix>.reference_ranking.tsv` records
  every candidate and its score. Congeneric samples skip the whole path — a same-genus reference is already the best
  obtainable — which keeps the extra NCBI calls off the large majority of a cohort. Controlled by
  `--enable_reference_reselection` (default true), `--n_reference_candidates`, `--reference_rank_read_subsample`.
  Degrades to the previous single-reference behaviour whenever a lookup fails or returns nothing usable.
- The Oatk fallback now routes on **length-inflated** MitoHiFi assemblies too, and the post-assembly routing rules are one
  predicate (`oatkFallbackReason`) rather than separate channels. A too-distant reference does not only truncate an
  assembly, it can inflate one: OG56 (*Vincentia punctata*, a *Fowleria vaiulae* reference at 76.2% identity) came back
  **circular with all 13 CDS** at 1.281x reference length, carrying a ~170 bp control-region unit in ~78.6 tandem copies.
  The gene count sees nothing wrong with that molecule, the empty-FASTA fallback never fires, the sample is
  `NON_CONGENERIC` rather than `CROSS_ORDER`, and `COLLAPSE_CONCATEMER` skips it because the anomaly is classified
  `control_region_repeat` and not `concatemer` — so it fell through every existing route and was published as-is. Length
  (from `MITOHIFI_CHECK_CIRCULARITY`'s `length_ratio`) is now checked as a second, independent measure of the same failure,
  at or above `--mitogenome_oatk_length_ratio_threshold` (default 1.15), and a `MISMATCH` relevance verdict routes on its
  own. `REFERENCE_RELEVANCE` is invoked inside the MitoHiFi subworkflow (aliased `REFERENCE_RELEVANCE_ROUTING`,
  unpublished) so the verdict exists early enough to route on.
  The length rule is **qualified** by that verdict and the gene rule deliberately is not: over-length alone is ambiguous,
  since OG848/OG852/OG853 all run 1.12-1.33x against congeneric references that `PASS`, which is real length heteroplasmy
  and not something a different assembler will fix — whereas OG2102 is `PASS` at 7 of 13 PCGs, so qualifying the gene rule
  would have un-routed the sample that branch was built for. Every rule demands positive evidence, so a missing or
  ungraded artefact leaves the assembly alone. Reasons are resolved before the channel mix, first match wins, because
  every routed sample is stamped with the same Oatk prefix and a sample arriving twice would launch two OATK tasks writing
  identical output names.
- `REFERENCE_RANK` no longer substitutes a reference it did not actually choose. All five of OG56's candidates scored
  `0.0000` with 0 mapped reads and the ranking still recorded `chosen=yes` against the first, so the assembly was built on
  an arbitrary member of an all-zero tie — while the module claimed re-selection "can only improve on the previous
  behaviour". A tie among zeros is the absence of a result, not a result. The three paths that reach a verdict without read
  evidence (no candidate recruits anything, no reads subsampled, a lone unmapped candidate) now mark every row `chosen=no`,
  record why in a new `status` column, and emit **empty** `chosen_reference` files; both assembly subworkflows test for
  size and keep `findMitoReference`'s original pick. This affected 8 samples in `mitogenomes-missing-audit-6`
  (OG56, OG81, OG750, OG801, OG2093, OG2202, OG2021, OG810). Deliberately **not** a routing signal: 7 of those 8 assembled
  acceptably, so an all-zero ranking does not predict a bad assembly — only the post-assembly gate above decides that.
- Gene-incomplete MitoHiFi assemblies now route to the reference-free Oatk fallback, not just empty ones. A reference too
  distant for read recruitment does not always yield *nothing*: it can drop the most divergent gene blocks and return a
  clean-looking but gene-incomplete molecule. OG2102 (*Rouleina attrita*, non-congeneric *Alepocephalus* reference) came
  back as a plausible 15.6 kb contig carrying 7 of 13 protein-coding genes, so the empty-FASTA fallback never fired. The
  PCG count is now read from MitoHiFi's own `final_mitogenome.gb` — available inside the assembly subworkflow, so no
  dataflow cycle — and anything below `--mitogenome_summary_expected_pcg_count` is sent to Oatk with
  `fallback_reason=gene_incomplete_mitohifi`. The MitoHiFi assembly stays published alongside the Oatk attempt so curation
  can compare rather than lose information.
- `family` and `order` are now resolved from the OceanOmics species table (`sp.family` / `sp.ordr`) and emitted onto the
  samplesheet (optional in `assets/schema_input.json`, so existing samplesheets stay valid). Without them every
  non-congeneric reference collapsed to the flat `NON_CONGENERIC` tier and the `CROSS_ORDER` route to reference-free
  assembly was unreachable dead code. `REFERENCE_DIVERGENCE` now also runs on the GetOrganelle reseed path, which
  previously produced no divergence flag at all.
- Post-assembly circularity re-check for MitoHiFi (`MITOHIFI_CHECK_CIRCULARITY` + `bin/check_circularity.py`). MitoHiFi's
  terminal-overlap test yields false negatives on hifiasm assemblies that are genuinely circular (the closed unitig loses
  its self-overlap once MitoHiFi rotates/trims it), which then mislabels the record in the SQL db as a "scaffold" and trips
  the assembly summary's `not_circularised` manual-review reason. For every sample MitoHiFi flags non-circular, the module
  remaps the already-mapped HiFi reads to a doubled reference and counts reads bridging the linearisation point
  (`--min-spanning-reads`/`--min-overhang`), and reads the hifiasm `c`/`l` contig flag; if either signal is positive it
  rewrites `was_circular` to `True` in the contig-stats table (a drop-in replacement consumed unchanged by the SQL upload
  and assembly summary). The original MitoHiFi verdict and all evidence are preserved in a published
  `<prefix>.circularity_check.tsv`. The same module also runs a length / tandem-repeat assessment on every assembly:
  it compares the assembly length to the related reference (parsed from the contig-stats header), self-aligns the
  assembly to locate any tandem-repeat array, and uses the final GenBank gene coordinates to decide whether that array
  sits in the control region. Over-length assemblies are classified as `concatemer` (tandem genome duplication — collapse
  to a monomer), `control_region_repeat` (D-loop VNTR / duplicated CR — review for length heteroplasmy) or `unresolved`
  (review for partial duplication / NUMT), with the redundant span stated as `suggested_trim_region` in the evidence file.
  This is a curation flag for manual review, never an automatic edit. The anomaly is also folded into the assembly
  summary's `manual_review_reason` (the `<prefix>.circularity_check.tsv` is read as a per-run sidecar), and the QC gate
  (`EVALUATE_QC_CONDITIONS`) blocks any sample with an anomaly from progressing to the submission-prep QC subworkflow
  (`proceed_qc = false`). Samples without a check (precomputed) get a `no anomaly` placeholder
  (`assets/placeholders/empty_circularity_check.tsv`) so they are never dropped or blocked on this condition.
- GetOrganelle circularity re-test + anomaly screen (`GETORGANELLE_CHECK` + `bin/check_getorganelle.py`). GetOrganelle
  reports a single non-circularised scaffold for assemblies it cannot formally close, but these are frequently complete
  circles linearised at a different origin (per the GetOrganelle docs). The check BLASTs the scaffold against the related
  reference (`relabel/<prefix>.reference.gb`): full reference coverage (default ≥95 %, single record, length within 0.9–1.15×)
  means a complete circle, and `meta.circular` is corrected to `true` so MITOS2 annotates it circular and the GenBank QC gate
  treats it as submission-ready; a contiguous uncovered chunk means a genuine gap and the verdict is left non-circular. The
  same module screens length / tandem-repeat anomalies (concatemer / tandem_repeat / unresolved) like the HiFi check. The
  corrected verdict + evidence flow to the assembly summary (`circularised` override, anomaly in `manual_review_reason`) and
  the QC gate (anomaly block; circular condition via `meta.circular`). Samples with no findMitoReference get an empty
  placeholder (`assets/placeholders/NO_REFERENCE.gb`) and are recorded as `no_reference` rather than dropped.
- GetOrganelle reseed now builds a custom gene (label) database from the reseed reference and passes it via `--genes`,
  improving recovery of divergent mitogenomes. New `GETORGANELLE_GENEDB` module + `bin/extract_getorganelle_genedb.py`,
  gated by `--getorganelle_genedb_min_genes` (default `10`).
- Per-sample mitochondrial genetic code derived from taxonomic `class` in samplesheet preparation
  (`meta.genetic_code`): Cnidaria → 4, echinoderms/flatworms → 9, other invertebrates → 4, vertebrates → the
  `--translation_table` default. Consumed by MITOS2, TRANSLATE_GENES, and MitoHiFi.
- Anthozoan annotation QC gate + reference-based fixer. MITOS2 annotates coral PCGs + 12S correctly but routinely
  drops the divergent 16S rRNA and one exon of the group-I-intron-split nad5. A new clade-aware gate
  (`ANNOTATION_QC_GATE` + `bin/annotation_qc_gate.py`) flags each invert annotation as FIX/PASS (deficient when 16S
  is missing or nad5 is truncated); only FIX corals are routed to `CORAL_ANNOTATION_FIX` (+ `bin/coral_fix_bed.py`),
  which BLAST-transfers the 16S and nad5 exons from a close coral reference into MITOS2's `result.bed` and re-runs
  `mitos_to_emma.py` (so the existing splice-join, translation, cds/proteins extraction and trnM re-origin are reused).
  Correctly annotated corals pass through MITOS2 untouched. The fixer is fail-safe (BLAST coverage/identity + a nad5
  ORF check guard every edit; a poor reference reproduces the original output). Reference is resolved per sample in
  priority order: the assembly stage's findMitoReference download → a fresh `MITOHIFI_FINDMITOREFERENCE` lookup →
  the bundled `assets/refdb/anthozoa/anthozoa_reference.gb`; the reference used is published into the sample's annotation dir.

#### `Fixed`

- A CDS that initiates on an alternative start codon is now declared as such even when Emma says nothing about
  it, so `SEQ_FEAT.StartCodon`/`SEQ_INST.BadProteinStart` stop quarantining otherwise clean assemblies.
  `bin/process_files.py` has injected `transl_except (pos:...,aa:Met)` for non-canonical starts since v1.1.0, but
  only for a CDS carrying Emma's `non-standard start codon` note. Emma writes that note inconsistently: it stayed
  silent on `OG663.ilmn.240313.getorg1770`, whose ATP6 initiates on `CTG`, so the note gate vetoed the fix before
  the sequence-level check could run and the sample failed the table2asn gate on those two errors alone. The note
  is no longer consulted as a trigger -- the codon read off the assembly is -- and where Emma did not write one the
  pipeline now writes `/note="non-standard start codon <CODON>"` itself so the flatfile records why the `aa:Met` is
  there.
  **The declaration is gated on the codon being a plausible initiator**, not merely on it being illegal under the
  sample's own table. NCBI translates the initiator as Met whatever it is, but only for codons its Starts row marks
  `M`, so an unrestricted rule would launder a mis-called CDS boundary or a single bad base call into a valid-looking
  record. `PLAUSIBLE_MITO_STARTS` is the union of the Starts rows across the mitochondrial translation tables, read
  out of the EMBOSS `EGC.*` data files shipped in the `seqret` container this pipeline already uses (tables 2, 3, 4,
  5, 9, 13, 14, 21; the later 24 and 33 add nothing new): `TTA, TTG, CTG, ATT, ATC, ATA, ATG, GTG`. `CTG` is a start
  in the mould/protozoan/coelenterate code, so OG663 clears; a start outside the union is left to fail validation and
  be looked at by hand, and is logged with the gene name and the offending codon rather than passing silently.
  Re-running the 20 assemblies of `batch-01` through the new code reproduces 19 feature tables byte for byte and
  changes only OG663's, whose regenerated `.val` is empty.

- `PREPARE_ENA_METADATA` no longer fails on every sample, so ENA candidate packages build again. The process
  read `ena_candidate_runs`, the last table of the in-repo ENA selection layer that migrations `010`-`012`
  retired; production never had it, so all 18 tasks of a run died with `relation "ena_candidate_runs" does not
  exist`. `errorStrategy = 'ignore'` kept the run alive, but the metadata channel is joined into
  `ch_ena_package_base`, so the join starved and `BUILD_ENA_CANDIDATE_PACKAGE` ran zero times. The read was dead
  weight regardless: nothing in the repo or the schema ever wrote a run accession, so the query could only ever
  return an empty list. `bin/prepare_ena_metadata.py` now emits `"run_accessions": []` as a constant, leaving the
  on-disk JSON schema and every `bin/ena_package.py` reader unchanged; `RUN_REF` is simply omitted from the
  manifest, which is the downstream submitter's to add since it owns the raw-read submissions to
  `PRJEB123419/420/421`. The `--run-accession` CLI escape hatch on `ena_package.py` stays.
  `sql/013_drop_ena_candidate_runs.sql` retires the table (no archive: it has never held a row), and the
  `selection_tables_dropped` audit in `bin/apply_ena_migrations.py` now covers it, which is why a schema the
  pipeline could not run against previously passed the post-migration audit clean.
- A sample with no recorded collection date no longer fails the table2asn gate. `bin/build_source_modifiers.py`
  filled a null `sample.date_collected` with the literal string `Unknown`, which table2asn rejects as
  `SEQ_DESCR.BadCollectionDate` ("Collection_date format is not in DD-Mmm-YYYY format"). The validation gate counts
  that as an `ERROR`, so one absent field quarantined an otherwise clean assembly -- `OG193.ilmn.240313.getorg1770`
  failed on this and nothing else. It is not a rare case: 626 of the 1747 sequenced `og_id`s (36%) have no date.
  The cell is now left empty, which omits the modifier, exactly as the `lat_lon` guard already does for a missing
  coordinate. `collection_date` is not a required source qualifier in the ENA flatfile and this pipeline validates
  with `-context sequence`, which registers no BioSample, so nothing downstream needs the field to be present.
  **The INSDC missing-value terms are deliberately not used here.** `missing`, `not collected` and `not provided`
  all pass table2asn and so look like the obvious fix, but they belong to the ENA *sample checklist* vocabulary
  (ERC000011), not to the flatfile qualifier: ENA's own `CollectionDateQualifierCheck` (sequencetools 2.33.2, as
  shipped in webin-cli 9.0.3) rejects all three, so adopting one would only move the failure from the first gate to
  the last. An empty cell is the only value that clears both. A `valid_collection_date()` guard now sits beside
  `valid_lat_lon()` and blanks anything that is not one of the three INSDC forms, logging the SeqID and the rejected
  value. It also catches **future** dates, which both validators refuse (`Collection_date is in the future` /
  `FutureDateException`): the `sample` table currently holds 60 day/month-transposed rows such as `2026-12-06`.
  None of those are sequenced yet, so nothing is blocked today, but they would have failed on assembly; correcting
  them is a database task, not a pipeline one.
- The ordered SQL migration chain can be replayed end to end again. `bin/apply_ena_migrations.py` reapplies every
  file on each run and relies on them being idempotent, but four of them referenced schema that a later migration
  had already removed, so a database part-way along the chain could not be brought forward at all -- the run aborted
  on `004` and left the schema where it was. Each is now guarded on the object it needs still existing, so a fresh
  database builds exactly as before while an advanced one skips the dead step:
  - `001` built `ena_validation_exact_result_idx` on `result_digest`, which `012` drops. Skipped when the column is
    gone; nothing is lost, since `002` drops that index and replaces it with `ena_validation_attempts_key_idx`.
  - `004` created `ena_candidate_loci` with a foreign key onto `ena_locus_registry (locus_tag)`, and `006` drops that
    column (`ERROR: column "locus_tag" referenced in foreign key constraint does not exist`). Skipped once the column
    is gone; `010` drops the table a few migrations later regardless.
  - `006` then altered `ena_candidate_loci` unconditionally, which fails once `010` has dropped it.
  - `011` dropped three columns from `ena_validation_attempts` without first dropping `ena_validation_latest`. That
    view is `SELECT DISTINCT ON (assembly_prefix) *`, and Postgres expands the `*` at creation time, so it depends on
    every column and blocks the drop. It is now dropped and rebuilt around the change, as `012` already did. `011`
    had never been run against a live database, so this was latent.
  The integration test's migration list also omitted `010`, which meant `012` was exercised against an ordering that
  cannot occur: `ena_candidate_loci` survived and its foreign key blocked the `ena_candidate_packages` drop.

- Assembly-summary rows now describe the assembly they are named for. Four defects in
  `bin/mitogenome_assembly_summary.py` shared one root cause -- association by unanchored substring test -- and between
  them corrupted 84 of the 250 rows in the `mitogenomes-missing-audit-6` table:
  - **Per-run sidecars no longer manufacture assemblies.** `strip_known_suffix` falls back to `Path(name).stem` for any
    suffix it does not know, so `<prefix>.reference_ranking.tsv` and `<prefix>.reference_candidates_status.tsv` each
    became an assembly of their own: **38 rows, 15% of the table**, every one reported `failed` for want of a FASTA and
    inflating the failure count from 16 real failures to 54. Both suffixes are now listed, and `discover_assembler_runs`
    also rejects structurally via `is_assembly_run_prefix` -- if an earlier dot-field carries the assembler token
    (`mitohifi`/`hifiasm`, `getorg`, `oatk`) and the last one does not, the tail is a suffix, not an assembly -- so the
    next sidecar nobody has written yet cannot reintroduce the bug. The file stays in the pool for `files_for_run`; it
    just no longer spawns a row.
  - **Gene counts are matched exactly, not by substring.** `parse_annotation_stats` joined with `og_id not in prefix` /
    `code not in prefix`. `OG5` is a substring of `OG58`, `OG8` of `OG810` and `OG848`, `OG10` of `OG107`, and
    `getorg1770` of `getorg1770reseed` -- so **45 rows reported counts they never earned** (OG5's 37 genes / 13 PCGs
    propagated across the cohort) and 2 rows that did have their own annotation reported another sample's numbers:
    `OG810.hic.260605.getorg1770reseed_rgj` showed 37/13/`no` where its own EMMA result was 25/10 with 12 missing genes.
    Because the loop returned the first match over an unsorted `rglob`, which sample won was filesystem-order dependent.
    The join is now on the filename EMMA writes (`<mt_assembly_prefix>.annotation_stats.csv`), falling back to whole-key
    equality rebuilt from the row's own `og_id.tech.seq_date.code`; a degenerate row with blank fields matches nothing
    rather than everything. `files_for_run` and `bin/multiqc_per_sample.py`'s `belongs_to_assembly` are anchored the same
    way -- the latter was copying the collapsed monomer's annotation, depth and LCA tables into the *pre-collapse*
    assembly's per-sample report (12 occurrences of `_collapsed` in OG750's).
  - **`num_genes` and `num_cds` come from one annotation or from neither.** `parse_mitohifi_stats` populated `num_genes`
    from `contigs_stats`' `number_of_genes`, which is MitoHiFi's own reference-guided annotation -- explicitly *not* the
    pipeline's gene count (see `subworkflows/local/mitogenome_assembly/mitohifi/main.nf`). OG750 therefore reported 56
    genes counted over its un-collapsed 2.14x concatemer beside a `num_cds` from the collapsed monomer. Both now come
    from the same `*.annotation_stats.csv` row, and an assembly that never reached annotation reports neither: blank
    means *not annotated*, which `is_complete_core` and `apply_qc` already treat as "not evaluated".
  - **A collapsed concatemer is reported as two honest rows instead of one blended one.** `COLLAPSE_CONCATEMER` renames a
    genuine collapse to `<prefix>_collapsed.fasta` and every downstream stage forks on that basename, so the monomer
    already had its own publish dir, `mitogenome_data` row and remap depth -- but its FASTA was never staged into the
    summary, leaving that row with no `final_length_bp` and a `failed` verdict, while `apply_collapse_override` wrote the
    monomer's length onto the pre-collapse row. `COLLAPSE_CONCATEMER.out.fasta` and `.out.evidence` now reach
    `ch_assembly_summary_files` (filtered to genuine collapses -- a passthrough emits `<prefix>.fasta`, which would
    collide in the module's flat staging dir), `collapse_child_evidence` reads the monomer's circularity from the
    post-curation check that is written under the *parent's* prefix, and the override is replaced by
    `apply_collapse_provenance`, which marks the original with the new `superseded` status so it is not triaged twice.
    OG750 goes from four wrong rows to two: `...v323mitohifi` `superseded` at its real 32672 bp, and
    `...v323mitohifi_collapsed` `manual_review` at 15293 bp, circular, 32 genes / 11 PCGs, missing `TS2;TD;CO2;TK;ATP8`.
- Reference-relevance check no longer flags good assemblies. It was calibrated on coral data (same genus ~99.6% identity,
  same family ~96.9%, wrong family ~81.5%) and applied unchanged to fish, whose mtDNA evolves far faster: in the
  `mitogenomes-missing-audit-5` run it called 27 assemblies `reference_mismatch`, **10 of them against a same-genus
  reference**, and held 18 finished mitogenomes (37 genes, 13 PCGs, circular, in-range length) at `manual_review`. The
  `PASS` and `MISMATCH` identity distributions did not separate, they abutted — PASS bottomed out at 88.4%, MISMATCH
  topped out at 87.7%. Four independent fixes, calibrated against all 95 assemblies in that run:
  - **Coverage is normalised by reference length, not assembly length.** Assembly-normalised coverage conflated a bad
    reference with an inflated or fragmented assembly: OG778 scored 0.30 against its *own species'* reference at 100%
    identity purely because the assembly was fragmented, and the control-region-repeat samples (OG852/OG853) were
    punished for being longer than any reference could cover.
  - **dc-megablast instead of blastn's default megablast**, matching the fix already made in `bin/check_getorganelle.py`
    (commit `e0d41f0`) but never carried across: megablast's long exact seeds miss diverged cross-species HSPs, reporting
    0.35 coverage on OG56 where dc-megablast reports 0.76.
  - **A congeneric reference is never called a mismatch** (capped at the new `DIVERGENT` state): it is the best reference
    obtainable, so low identity there is biology, not a labelling error.
  - **The identity floor is taxon-aware** — 82% for vertebrates, 88% (the validated coral value) for invertebrates.
  The verdict is now `PASS` / `DIVERGENT` / `MISMATCH` / `UNKNOWN`, and `MISMATCH` requires *both* poor coverage and low
  identity, since good coverage at low identity is a distant relative and high identity over part of the molecule is a
  partial assembly — neither is the wrong reference. On the audit-5 cohort this yields 0 mismatches and 14 advisory
  `DIVERGENT` calls, while the synthetic wrong-reference fixture still calls `MISMATCH`.
- `reference_mismatch`, `reference_divergent` and `no_congeneric_reference` are advisory on a complete-core assembly.
  They describe the *reference*, not the assembly; a circular molecule of the expected length with all 13 PCGs and both
  rRNAs is finished whatever reference built it. When a poor reference really did damage an assembly, the damage still
  blocks — `missing_protein_coding_genes` and every structural flag are unaffected — so OG2102, OG1946, OG56, OG675,
  OG696, OG769, OG810, OG852 and OG853 all remain in `manual_review`.
- Open-nomenclature species names are normalised to the ENA-submittable `Genus sp.` form, fixing
  `ERROR: Organism is not Submittable` rejections at webin-cli validation. ENA/NCBI only recognise `Genus sp.`
  for an undescribed species; `Genus sp` (no period) and `Genus spp.` are not taxa, and the flatfile's
  `/organism=` is taken verbatim from the nominal species name (e.g. `OG1834` was rejected on `"Chaunax sp"`,
  where `Chaunax sp.` is taxId 3041296 and submittable). New shared helper `bin/species_name_utils.py`
  (`normalise_open_nomenclature`) is applied in `bin/species_validation.py` — the path that actually reaches
  the flatfile, via `lca_results.tsv` → `evaluate_qc_conditions.py` → `FORMAT_FILES --species` → `process_files.py`,
  and which also populates `lca_validation.validated_species_name` — and in `bin/create_samplesheet.py` so the
  emitted `nominal_species_id` column agrees. `modules/local/validated_species_query/main.nf` inlines the same
  rule (its heredoc can't import from `bin/`) so the qc-only rerun path doesn't reuse a stale unnormalised name.
  Normalisation runs after the species/genus/family matching in `query_species_info()`, leaving
  `reference_species_id` (the MitoHiFi `findMitoReference` query) unchanged. Two intended knock-on effects in
  the QC gate: the `Found_in_blast_YN` substring test gets stricter (`chaunax sp` previously also matched
  `Chaunax spinosus`-style names), and `Genus spp.` samples — which never matched, since NCBI writes `sp.` —
  now pass the gate and proceed to QC. Postgres `sample.nominal_species_id` is left untouched as source of
  truth. Names needing a judgement call (`Centrodraco sp 2`, ``Nesogobius sp. `groove cheek` ``,
  `Synodus macrops cf`) are deliberately passed through unchanged.
- MITOS2 no longer hardcoded to genetic code 5: removed the `ext.code = 5` override that ignored per-sample taxonomy,
  so invertebrate (e.g. coral) annotations use the correct code.
- `--translation_table` is now the vertebrate/default fallback rather than a global override across all samples.
- `ena_validation_attempts` no longer grows a new row on every pipeline rerun. It previously deduplicated only on an
  exact `result_digest` match, so any rerun that changed so much as an error count or the flatfile hash (which most
  do) appended another history row under the same `--ena_validation_attempt` token. `push_ena_validation_results.py`
  now upserts on `(assembly_prefix, ena_study, validation_attempt)`: a rerun overwrites the previous attempt
  (`attempt_count` increments) as long as it hadn't reached `submission_ready`; once a row is submission-ready it is
  frozen, and a later rerun under the same token is reported as `locked` rather than overwriting the recorded
  success. `sql/002_ena_validation_attempts_single_row_per_attempt.sql` migrates existing tables (dedup down to one
  row per key, preferring a submission-ready row, then add `attempt_count`).
- The assembly summary's `status` column now means the same thing for every assembler, so it can be sorted,
  filtered and counted across a mixed cohort. MitoHiFi and Oatk rows already carried a computed QC verdict
  (`complete` / `manual_review` / `failed`), but GetOrganelle rows started from the tool's own log verdict and
  could terminate on `circular`, a value no other assembler could produce. All three now share one
  `finalise_status()` helper and one three-value vocabulary; topology stays in the `circularised` column
  rather than being encoded twice.
- GetOrganelle rows no longer reach `complete` on weaker evidence than the other assemblers. GetOrganelle
  reports either "circular genome" or "N scaffold(s)", and the scaffold form matched none of the parser's
  branches: `circularised` was left blank, so `not_circularised` never fired, nothing blocked, and an
  `unknown -> complete` promotion labelled the row `complete` despite the assembly never having been shown to
  be circular (8 rows in the `mitogenomes-missing-audit-5` cohort). `getorganelle_status_from_log` is
  replaced by `getorganelle_evidence_from_log`, which records `circularised=false` for any non-circular
  verdict; a `getorg_check.tsv` whose `final_verdict_circular` confirms circularity still wins over the log.
  This also removes a latent bug where `elif "complete" in status_text` was tested before
  `elif "incomplete" in status_text`: `"complete"` is a substring of `"incomplete"`, so a genuinely
  incomplete result would have been reported as complete.

#### `Dependencies`

#### `Deprecated`

- Locus-tag allocation is no longer this pipeline's job. Tags are assigned and injected by a separate
  downstream submission pipeline, so every file this one emits now carries **no** `/locus_tag` on any
  feature: the `.tbl`, the `.gbf`, the `.embl.gz` and the collaborator `.gff`. The qualifier is absent
  rather than blank because an empty `/locus_tag=""` fails both `table2asn` and Webin, so "blank" has
  exactly one submittable spelling. Removed: `ALLOCATE_ENA_LOCUS_TAGS` and `bin/allocate_ena_locus_tags.py`,
  the `<full_seqid>.locus_tag_mapping.tsv` artefact, the locus-tag rendering and GFF tagging in
  `bin/ena_package.py`, the mapping-to-`ena_candidate_loci` round trip in `bin/select_ena_submission.py`,
  and the `--ena_locus_prefix_{hifi,hic,ilmn}` params. `subworkflows/local/utils_ena_targets` still
  resolves the per-technology ENA child study, just not a prefix. Deleting the allocator loses nothing
  else: its `ena_specimen_accessions` upsert is duplicated in `select_ena_submission.py`, and its
  OG-numeric collision guard is already enforced by the `og_numeric NOT NULL UNIQUE` and
  `og_numeric = substring(og_id FROM 3)` constraints in `sql/004`. Two consequences worth knowing:
  `table2asn` now reports `FATAL: NO_LOCUS_TAGS` on every record, which stays advisory in
  `bin/parse_table2asn_validation.py` and must not be promoted back to fatal or it quarantines the whole
  run; and dropping the mapping from the package changes `package_digest` for every candidate.
  `sql/010_drop_ena_locus_tables.sql` retires `ena_locus_registry` and `ena_candidate_loci`, copying both
  into `*_archive` tables first, since the serials behind already-published tags cannot be rebuilt from
  the flat files. ENA-side locus-tag prefixes (`OGMTHIFI`, `OGMTHIC`, `OGMTILMN`) remain registered
  against the three child studies and are recorded in `assets/ena/accessions.tsv` for the downstream
  pipeline to use.

- `submission_ready` now means what this pipeline can actually attest, and the ENA selection layer is gone.
  The flag was hard-coded `false` by `bin/collate_ena_validation.py` and could only be earned later by a
  production Webin validation on an already-`SELECTED` package, which the pipeline never ran: every row ever
  written said `false`, including rows whose flatfile had passed every gate. It is now `true` exactly when the
  flatfile cleared this pipeline's last gate (`ena-webin-cli -context sequence`), which is the point at which a
  candidate is ready to hand to the submission pipeline. Selection and submission belong to that pipeline, so
  `ena_selection.nf` is removed along with `SELECT_ENA_SUBMISSION`, `WEBIN_VALIDATE_GENOME`,
  `RECORD_ENA_PACKAGE_VALIDATION`, `bin/select_ena_submission.py`, `bin/record_ena_package_validation.py`, and
  the `--ena_selection_mode` / `--ena_decision_file` / `--ena_selected_by` / `--ena_package_metadata` /
  `--ena_validate_webin_production` params. `sql/012_drop_ena_selection_layer.sql` drops
  `ena_candidate_packages`, `ena_submission_selections` and the `ena_submission_queue` view, and restores
  sql/001's original `CHECK (NOT submission_ready OR webin_status = 'PASS')`. The drop is irreversible from the
  repo: dump those relations first if their contents matter. Existing rows are deliberately not backfilled and
  are corrected on their next run.
- The ENA validation record drops from 46 columns to 27, ending at `submission_ready`. Removed: `package_status`,
  `webin_production_status` and the five `webin_*` production fields, `overall_status`, `package_digest`,
  `flatfile_name`/`_sha256`/`_size`, `manifest_name`/`_sha256`/`_size`, `workflow_run_name`,
  `workflow_session_id`, `pipeline_revision`, and `result_digest`. The production and package columns described
  work this pipeline does not do and always collated as `NOT_RUN`; `overall_status` collapses into
  `submission_ready` now that the two cannot disagree; the checksums and run provenance duplicate what the
  package's own `package_metadata.json` and the Nextflow run report already record. `push_ena_validation_results.py`
  also loses its cross-table freeze, which consulted the now-dropped `ena_submission_selections`: a rerun always
  overwrites its `(assembly_prefix, ena_study, validation_attempt)` row, and the `locked` outcome is gone.
  The two stale stub headers in `ENA_VALIDATION_RESULT` and `ENA_VALIDATION_SUMMARY`, which had drifted from the
  real output and would have failed the uploader's column check on any `-stub-run`, are regenerated from the
  column lists they mirror.
