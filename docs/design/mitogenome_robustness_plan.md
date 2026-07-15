# Mitogenome pipeline robustness plan

Status: in progress
Owner: Tyler Peirce
Source run analysed: `/scratch/pawsey1348/tpeirce/mitogenomes-missing-audit-3`

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
