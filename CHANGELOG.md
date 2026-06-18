# nf-core/oceangenomesmitogenomes: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/oceangenomesmitogenomes, created with the [nf-core](https://nf-co.re/) template.

### `Added`

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

- MITOS2 no longer hardcoded to genetic code 5: removed the `ext.code = 5` override that ignored per-sample taxonomy,
  so invertebrate (e.g. coral) annotations use the correct code.
- `--translation_table` is now the vertebrate/default fallback rather than a global override across all samples.

### `Dependencies`

### `Deprecated`
