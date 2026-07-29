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
