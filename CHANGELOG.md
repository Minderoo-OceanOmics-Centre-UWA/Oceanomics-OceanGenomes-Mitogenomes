# nf-core/oceangenomesmitogenomes: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/oceangenomesmitogenomes, created with the [nf-core](https://nf-co.re/) template.

### `Added`

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
  (`assets/empty_circularity_check.tsv`) so they are never dropped or blocked on this condition.
- GetOrganelle circularity re-test + anomaly screen (`GETORGANELLE_CHECK` + `bin/check_getorganelle.py`). GetOrganelle
  reports a single non-circularised scaffold for assemblies it cannot formally close, but these are frequently complete
  circles linearised at a different origin (per the GetOrganelle docs). The check BLASTs the scaffold against the related
  reference (`relabel/<prefix>.reference.gb`): full reference coverage (default ≥95 %, single record, length within 0.9–1.15×)
  means a complete circle, and `meta.circular` is corrected to `true` so MITOS2 annotates it circular and the GenBank QC gate
  treats it as submission-ready; a contiguous uncovered chunk means a genuine gap and the verdict is left non-circular. The
  same module screens length / tandem-repeat anomalies (concatemer / tandem_repeat / unresolved) like the HiFi check. The
  corrected verdict + evidence flow to the assembly summary (`circularised` override, anomaly in `manual_review_reason`) and
  the QC gate (anomaly block; circular condition via `meta.circular`). Samples with no findMitoReference get an empty
  placeholder (`assets/NO_REFERENCE.gb`) and are recorded as `no_reference` rather than dropped.
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
  the bundled `assets/anthozoa_reference.gb`; the reference used is published into the sample's annotation dir.

### `Fixed`

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

### `Dependencies`

### `Deprecated`

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
