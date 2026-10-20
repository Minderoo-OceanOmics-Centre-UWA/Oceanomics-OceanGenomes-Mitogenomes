# OceanOmics OceanGenomes Mitogenomes: Output

## Introduction

This page describes the files the OceanOmics OceanGenomes Mitogenomes pipeline publishes. All paths
are relative to the results directory (`--outdir`). Where a file feeds the MultiQC report, the
section it appears in is named so you can cross-reference the tables with the dashboard. The
standalone entrypoints `qc_only_from_annotations.nf` and `ena.nf` are covered at the end.

## Naming

Every assembly gets an **assembly prefix**, `<sample>.<sequencing_type>.<date>.<assembler><version>`,
for example `OG3092.ilmn.260925.getorg1770` (GetOrganelle 1.7.7.0), `OG910.hifi.241127.v3mitohifi`
(MitoHiFi 3) or `<sample>.hifi.<date>.v10oatk` (Oatk 1.0). Curated variants of an assembly get their
own prefix by suffix:

| Suffix | Meaning |
|---|---|
| `reseed` | GetOrganelle second pass seeded from a related reference ([details](#mitogenome-assembly-getorganelle)) |
| `_concat` | multi-contig assembly concatenated into one record for annotation ([details](#concatemer-collapse-and-fasta-sanitisation)) |
| `_rgj` | reference-guided join, off by default ([details](#concatemer-collapse-and-fasta-sanitisation)) |
| `_collapsed` | head-to-tail concatemer collapsed to one monomer ([details](#concatemer-collapse-and-fasta-sanitisation)) |

Each variant is a separate directory under `mitogenomes/<sample>/` with its own annotation, LCA and
QC results.

Annotation outputs add the annotator and its version, e.g. `.emma102` (EMMA 1.0.2) or `.mitos2110`
(MITOS 2.1.10). The prefix plus annotation version is the **`full_seqid`**
(`OG3092.ilmn.260925.getorg1770.emma102`). Everything from GenBank packaging onwards is named on it,
so two annotations of one assembly never overwrite each other.

## Results directory layout

```text
<outdir>/
├── samplesheet/                 samplesheet.csv, taxonomy_resolution.tsv
├── prepare/                     lca_cache_manifest.json
├── mitogenomes/
│   └── <sample>/
│       ├── fastp/               Hi-C read trimming reports (Hi-C only)
│       └── <assembly_prefix>/
│           ├── mtdna/           assembly, assembler logs, circularity and depth evidence
│           ├── annotation/      EMMA / MITOS2 annotation bundle and annotation QC
│           ├── annotation_rebuild/   only with --blast_precomputed_annotations
│           ├── lca/             BLAST hits and LCA calls per marker
│           ├── genbank/         table2asn inputs and outputs
│           ├── ena/             ENA metadata, candidate package and validation
│           └── <assembly_prefix>_multiqc_report.html
├── species_validation/          per-assembly species validation tables
├── qc/                          per-assembly QC summary
├── sql_uploaded_data/           database upload logs (only with --sql_config)
├── ena/                         run-level ENA readiness tables
├── multiqc/                     cohort MultiQC report and assembly summary
├── getorganelle_checkpoints/    GetOrganelle resume checkpoints
├── held_samples.tsv             samples held at any stage, with the reason
├── run_completeness.txt         every samplesheet sample accounted for
└── pipeline_info/               Nextflow reports, params, software versions
```

Runs launched with `nextflow_run.sh` also get `pipeline_info/launcher/` and `backup_scripts/`
([details](#pipeline-information)).

## Samplesheet preparation

<details markdown="1">
<summary>Output files</summary>

- `samplesheet/`
  - `samplesheet.csv`: the samplesheet the run used. In `--input_dir` mode this is the enriched
    sheet built by `CREATE_SAMPLESHEET_ENRICHED`, with columns `sample`, `sequencing_type`,
    `single_end`, `original_id`, `completion_date`, `date`, `assembly_prefix`,
    `nominal_species_id`, `reference_species_id`, `class`, `family`, `order`, `invertebrates`,
    `genetic_code`, `fastq_1` and `fastq_2`.
  - `taxonomy_resolution.tsv`: where each sample's lineage came from (the OceanOmics `species`
    table or the NCBI taxdump).
- `prepare/lca_cache_manifest.json`: checksum manifest for the shared FishBase and taxdump cache
  that every LCA task reads (see `--taxonkit_db_dir` in [usage.md](usage.md#external-resource-files)).

</details>

The samplesheet helper scans the staged FASTQ files, infers read pairing and sequencing technology,
and queries the database (`--sql_config`) for collection metadata and taxonomy. The taxonomic `class`
sets the per-sample genetic code and the annotation path. If the taxonomy check aborts the run, fix
`samplesheet/samplesheet.csv` by hand and rerun with `--input`.

## Mitogenome assembly (GetOrganelle)

Illumina and Hi-C samples are assembled with GetOrganelle. Raw Hi-C reads are trimmed with fastp
first (`--skip_hic_fastp` turns this off).

<details markdown="1">
<summary>Output files</summary>

- `mitogenomes/<sample>/fastp/`: `*.fastp.json`, `*.fastp.html`, `*.fastp.log` for Hi-C samples.
  Trimmed reads are not published.
- `mitogenomes/<sample>/<assembly_prefix>/mtdna/`
  - `<assembly_prefix>.fasta`: the assembly, FASTA header rewritten to the prefix.
  - `<assembly_prefix>.get_org.log.txt`: GetOrganelle log (parsed by the assembly summary).
  - `*.selected_graph.gfa`, `*extended_K*.assembly_graph.fastg`, `*.extend-*.fastg`,
    `*.extend-*.csv`: assembly graphs and their label tables.
  - `<assembly_prefix>.graph_join.tsv`: how the selected graph paths were combined.
  - `<assembly_prefix>.getorg_check.tsv`: circularity evidence from `GETORGANELLE_CHECK`, including
    any concatemer or tandem-repeat anomaly.
  - `<assembly_prefix>.mito_depth.tsv`: uniform read depth ([details](#read-depth)).
- `mitogenomes/<sample>/<assembly_prefix>reseed/`: the reseed pass, when one ran. Same `mtdna/`
  files as above, named `<assembly_prefix>reseed.*`, plus:
  - `mtdna/reference_seed/`: what seeded it.
    - Vertebrates: the NCBI reference as `<accession>.fasta` / `<accession>.gb`, the lookup result
      `<assembly_prefix>.findmitoreference_status.tsv`, and the custom gene database
      `<assembly_prefix>reseed.genedb.fasta` with its `.genedb_status.tsv`.
    - Invertebrates: `<assembly_prefix>.seed.fasta`, `<assembly_prefix>.reference.gb` and
      `<assembly_prefix>.genedb.fasta` from the curated group database, with the selection record
      `.reference_select.txt` (or `.fallback_select.txt` when the first pass was empty).
- `getorganelle_checkpoints/<assembly_prefix>_getorganelle/`: GetOrganelle working state kept so a
  `-resume` does not restart a long assembly from scratch.

</details>

GetOrganelle runs with `-R 20 -w 95 --continue` by default (override with
`--getorganelle_fromreads_args`).

**Automatic reseed.** When a first pass fails or comes back fragmented or non-circular because the
generic `animal_mt` seed was too divergent, the pipeline reruns GetOrganelle with a related
mitogenome as the seed (`-s`) and a gene database built from that reference as the labels (`--genes`).
Vertebrates take the reference from `MITOHIFI_FINDMITOREFERENCE`. Invertebrates take it from the
curated per-phylum database in `assets/refdb/<group>/`, chosen by BLASTing the first-pass assembly
against it, or by ranking candidates against a read subsample when the first pass produced nothing.
The vertebrate gene database is built only when the reference yields at least
`--getorganelle_genedb_min_genes` genes (default `10`); otherwise the sample keeps its first pass.
The better of the two assemblies is carried forward, so no sample is dropped. `--skip_getorganelle_reseed`
turns the reseed off.

## Mitogenome assembly (MitoHiFi)

<details markdown="1">
<summary>Output files</summary>

- `mitogenomes/<sample>/<assembly_prefix>/mtdna/`
  - `<assembly_prefix>.fasta`: circularised HiFi assembly, header normalised to the prefix.
  - `<assembly_prefix>.gb`: MitoHiFi's own GenBank annotation (not the one the pipeline uses).
  - `<assembly_prefix>.contigs_stats.tsv`: per-contig metrics (`MitoHiFi Contig Stats` in MultiQC).
  - `<assembly_prefix>.contigs_stats.with_coverage.tsv`: the same table with MitoHiFi's legacy
    `avg_coverage` / `coverage_cv` columns appended, corrected by the circularity check.
  - `<assembly_prefix>.coverage.tsv`: legacy MitoHiFi depth over reference-recruited reads only.
  - `<assembly_prefix>.circularity_check.tsv`: circularity evidence from `MITOHIFI_CHECK_CIRCULARITY`.
  - `<assembly_prefix>.mito_depth.tsv`: uniform read depth ([details](#read-depth)).
  - `<assembly_prefix>.reference_candidates_status.tsv`, `<assembly_prefix>.reference_ranking.tsv`:
    which NCBI references were considered and how they ranked.
  - `<assembly_prefix>.reference_divergence.txt`: the pre-assembly taxonomy guard (`CONGENERIC`,
    `NON_CONGENERIC`, `CROSS_ORDER`, ...).
  - `<assembly_prefix>.hifiasm.log`, `<assembly_prefix>.log`: hifiasm and MitoHiFi logs.
  - `all_potential_contigs.fa`, `shared_genes.tsv`, PNG plots, `contigs_circularization/`,
    `contigs_filtering/`, `MitoReference/`: MitoHiFi diagnostics. The read-mapping directories
    (`reads_mapping_and_assembly/`, `coverage_mapping/`) are not published.
- `assembly/`: a placeholder FASTA and `<prefix>.assembly_status.tsv` for any HiFi sample that
  produced no assembly at all, so it is still reported rather than silently dropped.

</details>

The prefix encodes the MitoHiFi version, so downstream files keep the provenance of the HiFi run.

## Mitogenome assembly (Oatk fallback)

<details markdown="1">
<summary>Output files</summary>

- `mitogenomes/<sample>/<assembly_prefix>/mtdna/` (prefix ends `.v10oatk`)
  - `<assembly_prefix>.fasta`: the structure-solved Oatk contig under the prefix name, so it flows
    into annotation like any other assembly.
  - `<assembly_prefix>.mito.ctg.fasta`: Oatk's native mitochondrial contig, kept for provenance.
  - `<assembly_prefix>.gfa`: Oatk assembly graph.
  - `<assembly_prefix>.oatk.log`, `<assembly_prefix>.oatk_status.tsv`: run log and outcome.
  - `<assembly_prefix>.getorg_check.tsv`: circularity evidence read from the GFA self-link.

</details>

These files exist **only** when `--enable_oatk_fallback true` and the MitoHiFi path could not be
used for the sample: no usable NCBI reference, a reference lookup that still failed after its
retries, a cross-order reference, or an empty MitoHiFi assembly. Oatk assembles the HiFi reads de
novo and finds the mitogenome by profile-HMM search against `--oatk_mito_db`, so it needs no species
reference. It is reported in the assembly summary as `assembler=Oatk`, with no reference species and
no `no_congeneric_reference` flag. See [usage.md](usage.md#reference-free-oatk-fallback) for setup.

## Read depth

`<assembly_prefix>.mito_depth.tsv` (in `mtdna/`, written for every assembler) is the only depth
figure comparable across GetOrganelle, MitoHiFi and Oatk, and across Illumina, Hi-C and HiFi. It
carries `mean_depth`, `median_depth`, `depth_cv`, breadth at 1x/10x/20x, the mitochondrial read
fraction, and the filter thresholds applied. It appears in MultiQC as `Mitogenome Read Depth` and
populates `mitogenome_data.mean_depth`. See [Mitogenome read depth](usage.md#mitogenome-read-depth)
for how it is measured.

## Mitogenome assembly summary

<details markdown="1">
<summary>Output files</summary>

- `multiqc/mitogenome_assembly_summary_mqc.tsv`: one row per assembly, shown in MultiQC as
  `Mitogenome assembly summary`.

</details>

`MITOGENOME_ASSEMBLY_SUMMARY` reads the assembly outputs above plus the annotation statistics.
Missing files are tolerated: the field is left blank and a warning goes to the task log.

Columns: `assembly_prefix`, `sample_id`, `assembler`, `status`, `final_length_bp`, `circularised`,
`num_candidate_contigs`, `num_final_contigs`, `num_genes`, `num_cds`, `missing_genes`,
`frameshift_flag`, `mean_coverage`, `coverage_cv`, `reference_species`, `reference_accession`,
`numt_flag` and `manual_review_reason`.

`status` is the QC verdict and means the same thing for every assembler, so the column can be sorted,
filtered and counted across a mixed cohort. It takes exactly three values: `complete` (a final
assembly exists and no blocking QC reason survived), `manual_review` (a final assembly exists but a
blocking reason stands, listed in `manual_review_reason`) and `failed` (the assembler errored or
produced no final assembly). It says nothing about topology; read `circularised` for that.

**Blocking and advisory reasons.** An assembly that is circular, in the expected length range and
carries all 13 protein-coding genes and both rRNAs is treated as a finished mitogenome (the
"complete core"). On a complete-core assembly the soft reasons (`low_mean_coverage`,
`high_coverage_variability`, `no_congeneric_reference`, `reference_divergent`, `reference_mismatch`,
and a tRNA-only gene shortfall) are recorded but advisory, which is why some `complete` rows still
carry text in `manual_review_reason`. Everything else blocks.

Manual review reasons are semicolon-separated and can include `missing_final_fasta`,
`not_circularised`, `multiple_candidate_contigs`, `multiple_final_contigs`, `missing_genes`,
`missing_protein_coding_genes`, `frameshift_detected`, `low_mean_coverage`,
`high_coverage_variability`, `possible_numt`, `ambiguous_getorganelle_graph`, `failed_run`,
`length_anomaly`, `length_outside_expected_range`, `reference_mismatch`, `reference_divergent`,
`no_congeneric_reference` and `data_limited`.

- `data_limited` marks assemblies limited by the input data (low coverage plus fragmentation) rather
  than by a pipeline defect.
- `missing_protein_coding_genes` fires when the annotation-derived CDS count is below
  `--mitogenome_summary_expected_pcg_count`, catching collapses that tRNA counts mask in the total.
- `ambiguous_getorganelle_graph` means GetOrganelle resolved more than one path through the graph.
- The reference reasons come from the pre-assembly reference-divergence guard
  (`no_congeneric_reference`) and the post-assembly reference-relevance check
  (`reference_divergent` / `reference_mismatch`). They describe the reference, not the assembly;
  when a poor reference actually damages an assembly, the damage itself (missing PCGs, or a
  structural flag) still blocks.

`reference_divergent` means the reference matches the assembly but is a distant relative, so a closer
reference would seed and annotate better. `reference_mismatch` is the stronger claim that the
reference neither covers nor matches the assembly, i.e. the species label pointed the reference
lookup at the wrong taxon. The check measures coverage against the **reference** length, so a
concatemer, a control-region repeat or a fragmented assembly does not read as a bad reference; it
never calls a congeneric reference a mismatch; and its identity floor is taxon-aware (82% for
vertebrates, 88% for invertebrates), because anthozoan mtDNA evolves far more slowly than teleost
mtDNA.

QC thresholds are set with `--mitogenome_summary_min_mean_coverage`,
`--mitogenome_summary_max_coverage_cv`, `--mitogenome_summary_min_length`,
`--mitogenome_summary_max_length`, `--mitogenome_summary_expected_gene_count` and
`--mitogenome_summary_expected_pcg_count` (defaults in
[usage.md](usage.md#assembly-summary-qc-thresholds)).

The parser can also be run by hand on an existing result tree:

```bash
mitogenome_assembly_summary.py \
  --input results/mitogenomes \
  --output mitogenome_assembly_summary_mqc.tsv
```

## Concatemer collapse and FASTA sanitisation

Three steps sit between assembly and annotation.

**`COLLAPSE_CONCATEMER`** runs only on assemblies whose circularity evidence says
`anomaly_type=concatemer`. It re-verifies the duplication by self-alignment and rewrites the assembly
as a single monomer; anything it cannot confirm passes through unchanged, so it can only ever fix a
genuine duplication. In the original assembly's `mtdna/`:

- `<assembly_prefix>.concatemer_collapse.tsv`: the collapse decision (`collapsed` or passthrough).
- `<assembly_prefix>.post_curation_check.tsv`: circularity evidence for the result.

A genuine collapse becomes a new assembly, `<assembly_prefix>_collapsed`, with its own annotation,
database row and depth measurement. Its `mtdna/` holds `<assembly_prefix>_collapsed.fasta` plus a
mirror of the original assembly's `mtdna/` files and the two reports above, so the curated directory
carries the full provenance next to its annotation. The assembly summary reports it at its monomer
length rather than flagging it as over-length.

Control-region tandem repeats are not collapsed. They stay in review with the suggested trim in the
circularity evidence, because the repeat may be real.

**`SANITISE_FASTA`** then prepares the FASTA for annotation. A single-contig assembly passes through.
A multi-contig assembly is concatenated into one record and published as
`<assembly_prefix>_concat/mtdna/<assembly_prefix>_concat.fasta`, so a partial mitogenome can still be
annotated and given a species ID.

**`GETORGANELLE_JOIN`** (`_rgj`, reference-guided join) is an alternative to that plain
concatenation and is **off by default** (`--skip_getorganelle_join true`). Where a multi-scaffold
assembly has a real reference, it BLASTs the scaffolds against it and joins them in reference order
with up to 100 bp of `N` between them. Output: `<assembly_prefix>_rgj/mtdna/<assembly_prefix>_rgj.fasta`
and `<assembly_prefix>_rgj.scaffold_join.tsv`.

Both `_rgj` and `_concat` are artificial molecules: neither is a submittable assembly. The join is
off because it never beat the concatenation: across the 12 NOVA_260909_LA samples where both forms
were annotated, it tied 10 and lost 2, each loss a gene landing on an `N` spacer (OG2998 `CYTB`,
OG3058 `ND2`). An inferred scaffold order adds nothing for annotation or species ID, while the joins
cost real gene models. Enable it with `--skip_getorganelle_join false` if you need the ordering for
another purpose.

## Mitogenome annotation

<details markdown="1">
<summary>Output files</summary>

- `mitogenomes/<sample>/<assembly_prefix>/annotation/`
  - `<assembly_prefix>.<annotation>.fa|gff|tbl`: the annotation (e.g. `.emma102`, `.mitos2110`).
    EMMA also writes an `.svg` map.
  - `cds/<GENE>.<assembly_prefix>.<annotation>.fa`: per-gene nucleotide FASTA. `CO1`, `RNR1` (12S)
    and `RNR2` (16S) are the BLAST inputs.
  - `proteins/MT-<GENE>.<assembly_prefix>.<annotation>.fa`: protein translations.
  - `<assembly_prefix>.annotation_stats.csv`: gene counts, missing genes, gene order and the PCG ORF
    check (`EMMA Annotation Stats` in MultiQC). This is the annotation QC verdict.
  - `<assembly_prefix>.reference_relevance.txt`: the post-assembly reference-relevance check.
  - `mtdna_rescue/`: status of the vertebrate ND4L/ATP8 and tRNA rescue passes, when one ran.
  - MITOS2 (invertebrates) only:
    - `mitos_raw/`: MITOS2's native output (`result.bed`, `.gff`, `.faa`, `.fas`, `.geneorder`,
      `.mitos`, `.png`, `.seq`).
    - `mitos_fix/`: the invertebrate/coral annotation fixer's status and corrected BED.
    - `mitos_conversion_audit.tsv`: how MITOS2 output was mapped to the EMMA layout.
    - `<assembly_prefix>.reference.gb`, `.seed.fasta`, `.genedb.fasta`, `.reference_select.txt`:
      the annotation reference picked label-free from the curated group database, when one was
      selected.
- `mitogenomes/<sample>/<assembly_prefix>/versions_emma.yml` or `versions_mitos.yml`: annotator
  version.

</details>

Vertebrate samples are annotated with EMMA; invertebrate samples (`invertebrates=true`) with MITOS2,
whose output is reshaped into the same layout as EMMA's. Both use the per-sample mitochondrial
genetic code resolved from taxonomic class during samplesheet preparation (e.g. Cnidaria → 4).

For vertebrates, two rescue passes patch the bundle in place when EMMA misses genes: one recovers
ND4L/ATP8 by alignment, one recovers tRNAs with tRNAscan-SE. Each publishes its status to
`annotation/mtdna_rescue/`.

The annotation bundle is reused by BLAST/LCA, the QC gate and GenBank packaging, so keeping the
whole directory lets you inspect or hand-curate an annotation and rerun from it
([Skipping annotation](usage.md#skipping-annotation)).

With `--blast_precomputed_annotations`, the bundle's `cds/` and `proteins/` are rebuilt from the GFF
and FASTA, and `annotation_rebuild/<annotation>.bundle_rebuild.tsv` records, gene by gene, whether
the bundle's own files matched the GFF.

## BLAST filtering and LCA

<details markdown="1">
<summary>Output files</summary>

- `mitogenomes/<sample>/<assembly_prefix>/lca/` (one set per marker: `CO1`, `12s`, `16s`)
  - `blast.<gene>.<assembly_prefix>.<annotation>.filtered.tsv`: hits with query span > 200 bp and
    identity at or above the gate (97% for invertebrates, 98% for vertebrates), with a timestamp and
    gene identifier appended.
  - `blast.<gene>.<assembly_prefix>.<annotation>.near.tsv`: hits from 90% up to the gate, for
    inspection only; nothing downstream reads it. If that band is empty it holds the 10 best hits by
    bitscore at any identity, so it is empty only when BLAST found nothing.
  - `filtered_summary.<gene>.txt`: per-gene hit summary.
  - `lca.<gene>.<assembly_prefix>.<annotation>.tsv`: lowest common ancestor call.
  - `lca_raw.<gene>.<assembly_prefix>.<annotation>.tsv`: the top hits (up to 20) that passed the
    LCA's own gates, with resolved taxonomy.
  - `lca_short.<gene>.<assembly_prefix>.<annotation>.tsv`: condensed LCA table, when produced.

</details>

The curated database (`--curated_blast_db`) is searched for vertebrates and `--nt_blast_db` for
invertebrates. The invertebrate search is restricted to Metazoa (`-taxids 33208`) so that `core_nt`
stays affordable, but not to the sample's own phylum: a contaminated or mislabelled sample must still
be able to hit outside it.

The search itself is unfiltered (the full hit list stays in the task work directory); the identity
gate is applied afterwards, and the LCA uses the same value. The LCA's coverage gate is
`--coverage_mode overlap --min_coverage 80`: aligned length as a share of the shorter of query and
reference. A query-cover gate would reject exact matches to barcode-length references when the query
is a whole gene.

## Species validation and QC

<details markdown="1">
<summary>Output files</summary>

- `species_validation/` (one set per assembly prefix)
  - `lca_results.<assembly_prefix>.tsv`: the LCA decision for each marker
    (`Species Validation (LCA Summary)` in MultiQC).
  - `lca_combined.<assembly_prefix>.tsv`, `blast_combined.<assembly_prefix>.tsv`: the combined
    tables passed to the uploads.
  - `taxonomic_evidence.<assembly_prefix>.tsv`: the per-region evidence the decision used.
  - `validation_record.<assembly_prefix>.json`: the record written to `lca_validation`.
- `qc/<assembly_prefix>.qc_summary.tsv`: species chosen for QC, annotation pass/fail and missing
  genes (`QC Summary` in MultiQC).
- `held_samples.tsv`: every sample held at any stage, with the stage and reason.
- `run_completeness.txt`: counts of samplesheet samples, held samples and submission-ready samples,
  and whether every samplesheet sample is accounted for. A sample missing from both lists means a
  task failed silently.

</details>

`SPECIES_VALIDATION` checks the sample's `nominal_species_id` against the markers. Its evidence for a
species- or genus-level label is the set of taxa the LCA grouped for each region (`species_in_LCA`,
the hits within 2 identity points of the top hit), plus any synonym or in-house `[organism=...]` name
on the BLAST rows inside that same window. A label that appears in the filtered table but well below
a different species' best hit is not evidence. The label is matched at the rank it actually names,
so a genus-only or family-only ID is validated at that rank. `EVALUATE_QC_CONDITIONS` then decides
whether the assembly proceeds to GenBank packaging.

Without `--sql_config`, the annotation statistics, species validation, QC gate, QC summary, held
samples report and assembly summary all still run, so a run without a database still reports gene
counts, missing genes and a completeness verdict.

## Database uploads

<details markdown="1">
<summary>Output files</summary>

- `sql_uploaded_data/`
  - `<assembly_prefix>.mtdna.upload.txt`: assembly row (`mitogenome_data`).
  - `<assembly_prefix>.annotation.upload.txt`: annotation statistics.
  - `<assembly_prefix>.lca_blast.upload.txt`, `<assembly_prefix>.lca_raw.upload.txt`: BLAST/LCA
    rows (`blast_filtered_lca`, `lca`, `lca_raw_results`).
  - `<assembly_prefix>.species_validation.upload.txt`: `lca_validation` row.
  - `<full_seqid>.ena_validation.upload.txt`: `ena_validation_attempts` row.
  - `<full_seqid>.qc_validator.upload.txt`: the second-validator write ([below](#automatic-second-validator)).
  - `upload_results_summary_mqc.tsv`, `upload_results_appendix.txt`: run summary across all uploads
    (`SQL Upload Summary` in MultiQC).

</details>

Uploads need `--sql_config` and are turned off by `--skip_upload_results`. Each upload writes a
readable log, so connectivity or data problems can be audited after the run. The ENA validation and
second-validator writes are non-fatal: a failure is logged and the run continues.

## GenBank and ENA packaging

Runs when `--sql_config` and `--ena_study` are both set and `--skip_ena_submission_prep` is not.

<details markdown="1">
<summary>Output files</summary>

- `mitogenomes/<sample>/<assembly_prefix>/genbank/`
  - `processed/<full_seqid>.fa|gff|tbl|cmt`: FASTA, GFF, feature table and structured comment with
    submission-ready headers (`FORMAT_FILES`).
  - `<sample>.bankit_metadata.csv`, `<sample>.bankit_metadata_latlon_cleaned.csv`: sample metadata
    from the database (`BUILD_SOURCE_MODIFIERS`).
  - `<full_seqid>.src`: source modifier table for table2asn.
  - `genes/`, `proteins/`: extracted gene sequences (`<full_seqid>.genes.fa` and per-gene files)
    and their translations, each header tagged with the genetic code.
  - `<full_seqid>.sqn`, `.val`, `.dr`, `.gbf` (and `.stats` when table2asn writes one): the
    table2asn package plus validator and discrepancy reports. `ERROR`, `REJECT` and discrepancy
    `FATAL` findings quarantine only that assembly; warnings are reported and it continues.
    `FATAL: NO_LOCUS_TAGS` is expected on every record and treated as advisory.
- `mitogenomes/<sample>/<assembly_prefix>/ena/`
  - `metadata/<full_seqid>.ena_input_metadata.json`: the inputs the package was built from.
  - `package/`: the ENA genome-context candidate, everything named on `full_seqid`: the EMBL flat
    file (`.embl.gz`), chromosome list (`.chromosome_list.tsv.gz`), feature table (`.tbl`),
    `package_metadata.json`, the collaborator handover files (`.fa`, `.gff`, `.genes.fa`) and
    `checksums.sha256`. [ena_package.md](ena_package.md) documents the package in full.
  - `validation/flatfile/`: the webin-cli sequence-context format check: `.webin_manifest.txt`,
    `.webin_status.tsv`, `.webin_validate.log`, `webin_output/`, and `validated/<full_seqid>.embl.gz`
    when it passed. This is the pipeline's last ENA gate.
  - `validation/<full_seqid>.ena_validation_result.tsv`: the result of every submission gate,
    including explicit `NOT_RUN`, `NOT_REQUESTED` and `NOT_APPLICABLE` states, with the annotation in
    its own column. `submission_ready` is `true` when the flat file cleared every gate this pipeline
    runs; that is the signal the downstream submission pipeline consumes.
- `ena/`
  - `ena_validation_results_mqc.tsv`: run-level readiness table (`ENA submission readiness` in
    MultiQC, alongside the `Table2ASN submission gate`, `Table2ASN findings`,
    `ENA flat-file conversion` and `ENA Webin validation` sections).
  - `ena_run_summary.tsv`: one row per assembly with each gate's result and the validation attempt
    token.

</details>

Table2asn or Webin failures do not stop a batch: the affected assembly is reported and left out of
the submission-ready set. The package records no locus tags (they are allocated downstream) and no
`STUDY`, `SAMPLE` or `RUN_REF`, because the study, BioSample and read submissions are registered by
the downstream submission pipeline. `ENA_SUBMISSION_PREP` only reads from the database, so it runs
under `--skip_upload_results` too. Its `COVERAGE` comes from this run's depth measurement, and the
package records which source it used as `mean_depth_source`.

### Automatic second validator

A mitogenome needs two validators recorded in `lca_validation` before it can be submitted.
`SPECIES_VALIDATION` fills the first (`validator = 'nf-core'`) when the nominal species ID is
supported by the LCA window of at least one region. `PUSH_QC_VALIDATOR` fills the second: any
assembly reaching `submission_ready = true` (it cleared table2asn, the EMBL conversion and the
webin-cli format check) gets `validator_2 = 'QCd-nf-core'`.

- **It never overwrites.** The `UPDATE` only fires when `validator_2` is NULL or blank, and there is
  no `--force`. A human reviewer who already signed off stays recorded.
- **It never inserts.** With no `lca_validation` row there is no first validator either, so the
  missing row is logged instead.

With `--ena_webin_validate false`, `submission_ready` is always `false`, so no second validator is
recorded. The per-assembly outcome (`validator_2_set`, `preserved`, `not_ready`, `no_row`,
`failed`) appears in the `Validator 2` column of the MultiQC SQL upload summary.

`ena_validation_attempts` keeps one row per `(full_seqid, ena_study, validation_attempt)`; a rerun
under the same key overwrites it. See [ENA Webin validation](usage.md#ena-webin-validation) for the
attempt token.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: the cohort report.
  - `multiqc_data/`, `multiqc_plots/`: parsed data and static plots.
- `mitogenomes/<sample>/<assembly_prefix>/<assembly_prefix>_multiqc_report.html`: one filtered report
  per assembly, built from the same inputs. Turn these off with `--skip_per_sample_multiqc true`.

</details>

Custom sections, defined in `assets/multiqc_config.yml`: `Mitogenome assembly summary`,
`Mitogenome Read Depth`, `MitoHiFi Contig Stats`, `Species Validation (LCA Summary)`,
`EMMA Annotation Stats`, `QC Summary`, `SQL Upload Summary`, `ENA submission readiness`,
`Table2ASN submission gate`, `Table2ASN findings`, `ENA flat-file conversion`,
`ENA EMBL preflight` and `ENA Webin validation`. The collated software versions are included too.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`,
    `execution_trace_<timestamp>.txt`, `pipeline_dag_<timestamp>.html`: Nextflow reports, one set
    per launch.
  - `params_<timestamp>.json`: the parameters of each launch.
  - `nf_core_oceangenomesmitogenomes_software_mqc_versions.yml`: every tool version.
  - `memory_hints.json`: the memory tier that succeeded per sample for GetOrganelle and Oatk, so a
    `-resume` requests it directly instead of climbing from the lowest tier again.
  - `pipeline_report.html`, `pipeline_report.txt`: only when `--email` / `--email_on_fail` is set.
  - `launcher/` (runs started by `nextflow_run.sh`): a snapshot of the launch script, the repo commit
    it ran from, and any `OG_*` overrides in the environment.
- `backup_scripts/backup.sh` (runs started by `nextflow_run.sh`): the backup script with this run's
  directory baked in. `sbatch backup_scripts/backup.sh` copies each `ena/package/` to
  `s3://ocom.oceangenomes/analysed-data/mitogenomes/curated/<sample>/<assembly_prefix>/` and the
  whole `mitogenomes/` tree to Acacia.
- `quarantine_stale_route/`: only when a rerun finds an assembly directory built under different
  routing (a change in the sample's class, invertebrate flag or genetic code, recorded in each
  assembly directory's `route.json`). The old directory is moved here rather than mixed with
  the new results, and the run log names it.

</details>

## Standalone QC-only workflow (`qc_only_from_annotations.nf`)

This entrypoint skips assembly, annotation, BLAST/LCA and the assembly uploads, and runs the
[GenBank and ENA packaging](#genbank-and-ena-packaging) stage on existing annotation files. It
publishes the same `genbank/` and `ena/` trees per assembly, the run-level `ena/` tables, and,
unless `--skip_upload_results true`, the `*.ena_validation.upload.txt` and
`*.qc_validator.upload.txt` logs.

The species name comes from `lca_validation.validated_species_name`, matched on fields parsed from
each file name (`og_id.tech.seq_date.code.annotation`), and the genetic code from the taxonomic
class that same query returns. A file with no database match gets species `unknown`.

## Standalone ENA runner (`ena.nf`)

`ena.nf` converts existing table2asn `.gbf` files to EMBL and validates them, or validates existing
EMBL files. It publishes per assembly under `mitogenomes/<sample>/<assembly_prefix>/ena/validation/`:

- `flatfile/`: as in the main pipeline.
- `local/<full_seqid>.ena_preflight_status.tsv`, `local/<full_seqid>.ena_preflight_check.tsv`: gzip
  and EMBL structure checks on input flat files (`--ena_mode validate` only).
- `<full_seqid>.ena_validation_result.tsv`.

plus `ena/ena_run_summary.tsv`, `ena/ena_validation_results_mqc.tsv` and a standalone MultiQC report.
With `--sql_config` and uploads on, it writes the same `*.ena_validation.upload.txt` logs.
