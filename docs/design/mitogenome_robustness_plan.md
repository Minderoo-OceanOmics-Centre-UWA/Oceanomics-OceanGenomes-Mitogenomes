# Mitogenome pipeline robustness plan

Status: implemented (Phase 2.1 Oatk needs a validation run)
Owner: Tyler Peirce
Source run analysed: `/scratch/pawsey1348/tpeirce/mitogenomes-missing-audit-3`

## Implementation status

| Phase | Change | State | Verification |
|---|---|---|---|
| 1 | Blocking-vs-advisory split + ambiguous-graph rewrite | Done | 15 samples manual_review -> done on the audit data, 0 regressions; 16 unit tests |
| 3.1 | Concatemer auto-collapse (script + module + summary + wiring) | Implemented; mixed-batch validation pending | breakpoint inferred from dominant whole-genome self-alignment; only concatemer evidence enters the module; OG750 shifted-boundary regression included |
| 3.2 | CR-VNTR / tandem repeat kept as review with trim surfaced | Done | trim already in circularity evidence; collapse passes these through |
| 4 | data_limited tag (low coverage + fragmentation) | Done | tags OG637/OG765/OG829 shallow attempts, leaves OG810 (26x) alone; unit tests |
| 2.2 | Reference divergence CROSS_ORDER detection and routing | Implemented; mixed-batch validation pending | CROSS_ORDER, no-reference and exhausted lookup-error reads bypass MitoHiFi into OATK; disabled OATK emits a structured data-limited record |
| 2.1 | Oatk reference-free HiFi fallback | **Complete + validated**, gated OFF | real Oatk assembly, module run via Nextflow+Singularity, summary parser, full-pipeline stub all pass |

The pre-audit unit suite had 52 passing tests and 3 BLAST-dependent concatemer
tests skipped on the login environment. The repair adds an OG750 shifted-boundary
regression and channel/process integration coverage; the mixed real-data batch is
still the final acceptance gate.

### Nextflow validation (run 2026-07-16, Nextflow 25.04.6)
- `nextflow run main.nf -preview`: full DAG builds with COLLAPSE_CONCATEMER, OATK
  (fallback enabled) and REFERENCE_DIVERGENCE wired in; "completed successfully".
- Standalone `-stub-run` of COLLAPSE_CONCATEMER and OATK: both execute, emitting
  the expected `.fasta` / `.concatemer_collapse.tsv` / post-curation evidence /
  `.mito.ctg.fasta` / `.oatk.log` / OATK status TSV.
- Full-pipeline `-stub-run` (mitohifi path, Oatk fallback enabled): MITOHIFI_MITOHIFI,
  REFERENCE_DIVERGENCE, OATK and MITOGENOME_ASSEMBLY_SUMMARY all COMPLETED; the
  assembly-summary MultiQC table was produced. The only FAILED tasks were the two
  DB-download stubs (DOWNLOAD_BLAST_DB / DOWNLOAD_TAXONKIT_DB, exit 127 in stub),
  which the pipeline error-ignores and which are unrelated to these changes.
- Note: the repo's config requires Nextflow 25.x; the config parser in 26.04 rejects
  it (unrelated to these changes).

### Running the Oatk fallback (Phase 2.1)
Everything is wired; enable it in two steps.

1. Fetch the OatkDB profile-HMM for your clade (ray-finned fish shown), once:
   ```
   bash bin/download_oatk_db.sh actinopterygii_mito /scratch/$USER/oatk_db
   ```
   This pulls `<clade>_mito.fam` plus its `.h3f/.h3i/.h3m/.h3p` index files.

2. Run the pipeline with the fallback on:
   ```
   nextflow run main.nf <your usual opts> -profile singularity \
       --enable_oatk_fallback true \
       --oatk_mito_db /software/projects/pawsey0964/oatk_db/actinopterygii_mito.fam
   ```

Container: `quay.io/biocontainers/oatk:1.0--h577a1d6_1` (bundles oatk + syncasm +
nhmmscan/hmmer). Params (`enable_oatk_fallback`, `oatk_container`, `oatk_mito_db`,
`oatk_syncmer_size`, `oatk_syncmer_coverage`) default in `nextflow.config`'s `params`
block; the OATK resource request lives in `conf/modules.config` alongside every other
process. No dedicated `-c` config: the feature is gated by `enable_oatk_fallback`.

Behaviour: runs only on samples MitoHiFi failed to assemble; Oatk contigs are
named `<id>.<seqtype>.<date>.v10oatk.fasta`, annotated by the same EMMA/MITOS2
path, and reported in the assembly summary as `assembler=Oatk` (reference-free, so
no reference species / no_congeneric flag). Circularity is read from the Oatk GFA
self-link.

### Validation performed (2026-07-16)
- Container: oatk 1.0 + nhmmscan (HMMER 3.4) run; actinopterygii_mito.fam readable.
- Real assembly: simulated fish-mito HiFi reads -> Oatk -> 16,465 bp single circular
  contig (exact match to the source mitogenome).
- Module via Nextflow + Singularity (real container + DB): emits `<prefix>.fasta`
  (16,465 bp) + `<prefix>.gfa` (circular self-link). Caught + fixed a DB-staging bug
  (the `.fam` index files must be staged beside the `.fam`).
- Summary parser: Oatk run -> status=complete, circularised=true, 37/13 genes,
  reference-free; unit-tested (8 Oatk tests; 37 total green).
- Full-pipeline `-stub-run` (Oatk enabled): MITOHIFI -> OATK -> SANITISE -> EMMA ->
  SUMMARY all complete.

Commits: see branch `mitogenome-robustness` (5 commits, one per phase). The
in-progress ENA work in the tree was left untouched (never staged).

## Baseline

142 unique samples. Best result per sample across all assembly attempts:

| Best outcome | Samples |
|---|---|
| circular | 63 |
| complete | 50 |
| manual_review | 24 |
| failed | 5 |

113 / 142 (80%) auto-complete. The 29 that do not are the target.

### Confirmed non-causes
- **LCA / species-validation never gates a sample.** Zero manual-review reasons are species/`reference_mismatch`. Validation ran cleanly for every sample.
- **Annotation is adequate everywhere** except a handful of single-tRNA shortfalls (all 13 PCGs + 2 rRNAs present).

All real blockers are in **assembly-QC classification** and **MitoHiFi reference-based read recruitment**.

### Guiding principle (Phase 1)
An assembly that is circular, has all 13 PCGs + 2 rRNAs, and is in the expected length range is a finished mitogenome. Soft signals (coverage, no-congeneric reference, graph shape, one missing tRNA) should annotate it, not block it.

---

## Root-cause buckets

| Bucket | Cause | Representative samples |
|---|---|---|
| 1 | `ambiguous_getorganelle_graph` false positive on clean single-path circular assemblies (flag keys off GFA segment count, not path count) | OG100, OG12, OG1791, OG1849, OG56, OG61, OG658, OG698, OG793, OG804 |
| 2 | MitoHiFi recruits 0/too-few reads because `findMitoReference` returned a divergent (cross-genus / cross-order) reference; hifiasm gets nothing | OG62, OG109, OG1422, OG2093 (+ OG2102 partial collapse) |
| 3 | `no_congeneric_reference` forces review on otherwise-complete assemblies | OG1369, OG2104, OG2202, OG801 |
| 4 | `low_mean_coverage` forces review on complete circular 37/13 assemblies | OG1161, OG766, OG778, OG819, OG863 |
| 5 | Over-length repeat structures with the trim already computed but never applied | OG750 (concatemer), OG852, OG853, OG1946 |
| 6 | Single/double missing tRNA forces review | OG675, OG696, OG756, OG848 |
| 7 | Genuinely reads-limited (HiC low mito content / low-yield Illumina) — not a pipeline defect | OG765, OG637, OG829, OG810, OG864 |

---

## Plan

### Phase 1 — Classification fixes (`bin/mitogenome_assembly_summary.py`, Python only)

Introduces a **blocking vs advisory** split: all reasons stay in `manual_review_reason` for transparency, but advisory reasons no longer flip `status` to `manual_review` when the assembly passes a `complete_core` guard (circular + all PCGs + length in range + at most `TRNA_TOLERANCE` missing tRNAs).

- **1.1** Rework `getorganelle_graph_ambiguous()` to key off the count of resolved `path_sequence.fasta` files (>1 = genuinely ambiguous), not raw GFA `S`-segment count.
- **1.2** `no_congeneric_reference` advisory when `complete_core`.
- **1.3** `low_mean_coverage` advisory when `complete_core`.
- **1.4** `missing_genes` advisory when the shortfall is tRNA-only (all PCGs present, shortfall ≤ `TRNA_TOLERANCE`). `missing_protein_coding_genes` stays blocking.

Verification: re-run the summary generator over the audit run outputs and confirm the reclassification.

### Phase 2 — Reference-free HiFi assembly recovery

- **2.1** Add **Oatk** (reference-free, profile-HMM, HiFi-native) as a fallback branch fed by the existing empty/failed MitoHiFi FASTA branch. Requires a new module, container, and an actinopterygian mito OatkDB.
- **2.2** Harden the reference finder: reject picks outside the sample's order/family (e.g. OG1422 beachsalmon → flatfish); on rejection route to the Oatk fallback. Optionally rank candidates by sketch distance to a read subsample rather than taxonomy rank.

### Phase 3 — Auto-curation of over-length assemblies

- **3.1** Auto-collapse clean ~2x concatemers using the `suggested_trim_region` already computed by the circularity check, then re-annotate.
- **3.2** Control-region VNTR / tandem repeat: keep as review (possible real heteroplasmy) but surface the suggested monomer for one-click curation.

### Phase 4 — Data-limited labelling

Label reads-limited samples (HiC low mito content, low-yield Illumina) distinctly from pipeline failures so they are not counted as pipeline problems. No assembler change recovers these; route to a better library type where one exists.

---

## Projected outcome

| Phase | Recovered | Cumulative |
|---|---|---|
| Baseline | - | 113 / 142 (80%) |
| Phase 1 | ~15 | ~128 (90%) |
| Phase 2 | ~5 | ~133 (94%) |
| Phase 3 | ~1-4 | ~134-137 (~96%) |
| Phase 4 | 0 (relabel) | residual ~5 correctly classified |

## References
- Oatk (Genome Biology 2025): https://link.springer.com/article/10.1186/s13059-025-03676-6
- Oatk README: https://github.com/c-zhou/oatk
- MitoHiFi (2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10354987/
