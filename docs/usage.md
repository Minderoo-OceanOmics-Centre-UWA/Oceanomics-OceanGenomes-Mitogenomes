# OceanOmics OceanGenomes Mitogenomes: Usage

The OceanOmics OceanGenomes Mitogenomes workflow assembles, annotates, validates, and packages
mitochondrial genomes from PacBio HiFi, Illumina, and Hi-C sequencing runs. This page describes
pipeline-specific input requirements and the parameters you must provide to execute a full run on
OceanOmics data holdings.

> For general guidance on Nextflow execution, configuration profiles, and infrastructure tuning,
> refer to the nf-core documentation: <https://nf-co.re/docs/usage>

## Quick start

```bash
nextflow run Minderoo-OceanOmics-Centre-UWA/Oceanomics-OceanGenomes-Mitogenomes \
  -profile singularity \            # or docker / conda / podman / …
  --input oceanomics_samples.csv \   # or --input_dir "/data/OG*/**/*.fastq.gz"
  --outdir results \
  --organelle_type animal_mt \
  --curated_blast_db /path/to/OceanGenomes.CuratedNT.fasta \
  --nt_blast_db /path/to/nt/db/core_nt \
  --sql_config /path/to/postgres.cfg \
  --blast_db_dir /scratch/databases/blast \
  --taxonkit_db_dir /scratch/databases/taxonkit \
  --template_sbt /path/to/template.sbt
```

The command above runs the full assembly + annotation + validation + QC pipeline. If you already
have annotation files and only need QC packaging, use the standalone workflow below.

## Smoke test (`-profile test`)

```bash
nextflow run . -profile test,singularity -stub --outdir results_test
```

Runs two bundled invertebrate HiFi rows through samplesheet parsing and the assembly
wiring, in stub mode, with no external inputs. Annotation, ENA prep and the SQL
upload/QC modules are skipped, because each needs a resource this repository cannot
ship (`--taxonkit_db_dir`, `--nt_blast_db`, `--template_sbt`, Postgres credentials).
Use it to check that a clone and its container engine work, not to produce data.

## Standalone QC-only workflow (`qc_only_from_annotations.nf`)

```bash
nextflow run qc_only_from_annotations.nf \
  -profile singularity \
  --annotation_files "/path/to/mitogenomes/*/*/annotation/*.{fa,fasta,gff,tbl,gb}" \
  --sql_config /path/to/postgres.cfg \
  --template_sbt /absolute/path/to/template.sbt \
  --outdir qc_results
```

Required parameters for this standalone workflow:

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--annotation_files` | ✔ | Glob to precomputed annotation files (`.fa/.fasta/.gff/.tbl/.gb`) for each assembly prefix. |
| `--sql_config` | ✔ | PostgreSQL config used to fetch `validated_species_name` from `lca_validation` and to build source modifier tables. |
| `--template_sbt` | ✔ | `.sbt` template passed to `table2asn` via `GEN_FILES_TABLE2ASN`. |
| `--outdir` | ✔ | Output directory for published GenBank/QC artefacts. |
| `--translation_table` | Optional | Fallback mitochondrial genetic code (default `2`). The code is normally resolved per sample from the taxonomic class the SQL query returns; this value is used only for a class with no confirmed code, and the run warns for each such sample. |
| `--skip_upload_results` | Optional | Skip the ENA validation and second-validator writes (default `false`). |

Filename parsing and SQL matching rules:

- Input basenames must include at least five dot-delimited fields:
  `<og_id>.<sequencing_type>.<seq_date>.<code>.<annotation>[...]`.
- The workflow maps parsed fields to SQL keys as:
  `og_id = part[0]`, `tech = part[1]`, `seq_date = part[2]`, `code = part[3]`,
  `annotation = part[4]`.
- Query target: `lca_validation.validated_species_name`.
- Retrieved `validated_species_name` is passed into `FORMAT_FILES` and becomes the species value in
  processed FASTA/GFF headers used for submission files.
- If no row matches, species defaults to `unknown`.

Nextflow ≥ 24.10.x is required. Always combine the pipeline with an execution profile (`-profile`)
that matches your container/conda environment.

## Essential parameters (main pipeline)

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--input` | ✔ (or `--input_dir`) | Path to a CSV samplesheet that matches `assets/schema_input.json`. |
| `--input_dir` | ✔ (or `--input`) | Glob pointing at FASTQ files when you want the pipeline to build the samplesheet automatically. |
| `--outdir` | ✔ | Destination for all published results. |
| `--organelle_type` | ✔ | Organelle label passed to GetOrganelle (e.g. `animal_mt`). |
| `--curated_blast_db` | ✔ | NCBI-format BLAST database used for species validation. Provide an absolute path. |
| `--nt_blast_db` | Conditional | NCBI nt BLAST database PREFIX (e.g. `.../db/core_nt`) used when samples are marked `invertebrates=true`; searched with `-taxids 33208`. Must be nt-scale — the RefSeq-only `mito` database yields almost no invertebrate hits. |
| `--mitos_refdb` | Conditional | MITOS2 RefSeq reference data directory (parent of the version dir), downloaded from Zenodo. Required whenever invertebrate samples (annotated with MITOS2) are present. |
| `--mitos_refseq_ver` | Optional | MITOS2 reference version subdir name passed to `runmitos.py -r` (default `refseq89m`). |
| `--sql_config` | Conditional | INI file with `[postgres] dbname,user,password,host,port`. Required for `--input_dir`. Optional for `--input`: without it SQL uploads and ENA submission prep are off, and local QC still runs. |
| `--blast_db_dir` | ✔ | Directory used to cache the downloaded `taxdb.*` files; re-use between runs to avoid repeated downloads. |
| `--taxonkit_db_dir` | ✔ | Directory used to cache the NCBI taxdump for TaxonKit. |
| `--template_sbt` | ✔ | Submission template passed to `table2asn` when packaging GenBank artefacts. |
| `--ena_webin_validate` | Optional | Format-validate each EMBL flatfile with `ena-webin-cli -context sequence` (default `true`). This is the pipeline's last ENA gate; it needs no BioSample and never submits. |
| `--webin_validate_timeout_seconds` | Optional | Per-attempt wall-clock bound on the `webin-cli -validate` call (default `900`). A healthy call returns in seconds; the bound stops a hung ENA request from consuming the task's whole allocation. |
| `--webin_validate_max_attempts` | Optional | How many times `WEBIN_VALIDATE` calls `webin-cli` before giving up (default `3`). Only timeouts and infrastructure failures are retried; a validation failure is deterministic and stops immediately. |
| `--ena_study` | ENA validation | ENA study the run validates against, recorded in each package as `validation_study`. There is no default. It is not a submission target: the study an assembly is submitted into is the BioProject the downstream submission pipeline registers, so this value never reaches a submission manifest. |
| `--samplesheet_prefix` | Optional | Reserved for generated samplesheet naming in wrapper scripts. |
| `--getorganelle_genedb_min_genes` | Optional | Minimum genes a reference must yield to build the reseed custom gene database (default `10`). Below this, the sample keeps its first-pass GetOrganelle assembly instead of reseeding. |
| `--getorganelle_fromreads_args` | Optional | Override the default GetOrganelle from-reads arguments (default `-R 20 -w 95 --continue`). |
| `--enable_oatk_fallback` | Optional | Enable the reference-free Oatk HiFi fallback (default `false`). See [Reference-free Oatk fallback](#reference-free-oatk-fallback). |
| `--oatk_mito_db` | With Oatk fallback | Path to the OatkDB `<clade>_mito.fam` profile-HMM (index files staged beside it). Required when `--enable_oatk_fallback true`. |
| `--oatk_syncmer_size`, `--oatk_syncmer_coverage` | Optional | Syncmer size (`-k`, default `1001`) and coverage (`-c`, default `30`) passed to Oatk/syncasm. |
| `--skip_getorganelle_reseed` | Optional | Disable the automatic GetOrganelle reseed pass on failed/fragmented first passes (default `false`). |
| `--getorganelle_firstpass_group_genes` | Optional | Label the FIRST GetOrganelle pass with the sample's curated group gene database instead of only the stock `animal_mt` LabelDatabase (default `true`). The stock seed still recruits the reads; only labelling changes. |
| `--getorganelle_empty_first_pass_rescue` | Optional | Reseed a first pass that produced no contig at all from the sample's whole curated group database (default `true`). Without it an empty first pass cannot select a seed subset and can never be reseeded. |
| `--hinted_memory_base_gb` | Optional | Attempt-1 memory (GB) for `GETORGANELLE_FROMREADS`, `GETORGANELLE_RESEED` and `OATK` when no per-sample hint exists (default `64`). On Setonix the memory request is what is billed, so this dominates SU. |
| `--skip_hic_fastp` | Optional | Skip fastp trimming of raw Hi-C reads before GetOrganelle (default `false`). Override the fastp arguments with `--hic_fastp_args`. |
| `--force_db_overwrite` | Optional | Replace an existing `mitogenome_data` row on upload (default `true`). Set `false` for insert-only uploads, which keep a row that already holds a real result. Neither mode replaces a successful assembly row with `failed to assemble`. |
| `--lca_prune_history` | Optional | Delete `lca` / `lca_raw_results` rows this run did not reproduce (default `false`). Those tables are content-keyed upserts, so a rerun is correct without it and superseded results are kept as history. |
| `--translation_table` | Optional | Mitochondrial genetic code for vertebrate/unresolved samples (default `2`). Invertebrate codes are derived per sample from the taxonomic `class` via `assets/taxonomy/mito_genetic_codes.json` (Cnidaria, Ctenophora, Porifera → 4; Mollusca, Arthropoda, Annelida and most other protostomes → 5; Echinodermata, Platyhelminthes → 9; Tunicata → 13), so this does not force a single code across the run. |
| `--allow_unknown_taxonomy` | Optional | Continue with a warning when a sample's taxonomic `class` is `unknown` (default `false`, i.e. the run aborts). An unknown class is treated as vertebrate by every downstream chooser, so this is for deliberate one-offs only. |

Invertebrate reseeds take their seed and gene-label database from the curated per-phylum databases in
`assets/refdb/<group>/`, chosen from the sample's `class` (Cnidaria → anthozoa, Gastropoda → mollusca,
Malacostraca → arthropoda, and so on). An invertebrate class no database covers is **not** reseeded from another
phylum: it keeps its first-pass assembly, and the run log names it. See
[`assets/refdb/README.md`](../assets/refdb/README.md) for what each database contains and how to add a group.

The resolved per-sample code reaches every artefact that states one: the genome FASTA, the
`.tbl`, the `mgcode` passed to table2asn, the translation table, and the `[mgcode=]` tag on
each extracted gene, CDS and protein header. A sample whose `class` is marked
`invertebrates=true` but has no confirmed code aborts the run rather than falling back --
set the `genetic_code` column for it (commonly 5 for Mollusca/Annelida/Arthropoda/Porifera,
9 for Echinodermata/Platyhelminthes, 13 for Tunicata). Note that the invertebrate annotation
path is tuned for corals (cox1 origin rotation, coral reference DB), so review the annotation
of any non-cnidarian sample you push through it.

> **Assemblies published before this became per-sample.** Gene extraction used to hardcode
> `[mgcode=2]`, so any non-code-2 assembly QC'd before that fix has extracted genes, CDS and
> proteins that contradict its own genome record -- including inside the ENA candidate
> package. The genome FASTA, `.tbl` and table2asn output were always correct; only the
> extraction artefacts are wrong. Find them with:
>
> ```bash
> # an assembly whose genome record and gene records disagree
> grep -l 'mgcode=2' <outdir>/mitogenomes/*/*/genbank/genes/*.genes.fa
> ```
>
> then re-run those through `qc_only_from_annotations.nf`, which re-derives the code per
> sample from the SQL taxonomic class, so no `--translation_table` is needed for a mapped
> class:
>
> ```bash
> nextflow run qc_only_from_annotations.nf \
>   --annotation_files '<outdir>/mitogenomes/OG####/*/annotation/*.{fa,fasta,gff,tbl,gb}' \
>   --sql_config <cfg> --template_sbt <sbt> --outdir <outdir> \
>   --skip_upload_results true -profile setonix
> ```
>
> Run the first batch with `--skip_upload_results true`, confirm the headers agree, then
> re-run without it if the ENA validation attempts should be recorded. Any sample the run
> warns about with "no confirmed mitochondrial genetic code for class" fell back to
> `--translation_table` and needs its code passed explicitly.

`--input_dir` mode requires `--sql_config` because the enriched samplesheet generator queries
OceanOmics metadata. Pass `--taxonkit_db_dir` as well: the generator falls back to the NCBI
taxdump for samples the `species` table cannot resolve, and without it those samples fail the
taxonomy check. When running with `--input`, you can omit `--sql_config`: SQL uploads and ENA
submission prep are then skipped with a warning, and local QC still runs.

### Reference-free Oatk fallback

MitoHiFi recruits reads by mapping to the closest NCBI reference returned by `findMitoReference`. For lineages with
no close reference, it can recruit ~zero reads and produce nothing. The optional Oatk fallback (`--enable_oatk_fallback
true`) instead assembles the HiFi reads de novo and locates the mitogenome by profile-HMM search, so it needs no
species reference. It runs when no usable reference is found, a technical lookup still fails after three retries,
the selected reference is cross-order, or MitoHiFi emits no assembly. These outcomes remain distinct in the status
TSVs so a biological reference gap is not confused with an infrastructure failure.

When enabled, `--oatk_mito_db` is **required**: point it at the OatkDB `<clade>_mito.fam` profile-HMM for the sample's
clade, with its `.h3f/.h3i/.h3m/.h3p` index files staged in the same directory. Fetch a database once with the bundled
helper:

```bash
bash bin/download_oatk_db.sh actinopterygii_mito /scratch/$USER/oatk_db
# then: --enable_oatk_fallback true --oatk_mito_db /scratch/$USER/oatk_db/actinopterygii_mito.fam
```

Browse available clade databases at <https://github.com/c-zhou/OatkDB>. The Oatk container bundles `oatk`, `syncasm`,
and `nhmmscan`.

### ENA Webin validation

Every viable assembly/annotation version produces a self-contained ENA
genome-context candidate package during the sequencing run. The full OceanOmics
SeqID is retained in the EMBL entry, chromosome list, package metadata, and filenames.
The chromosome-list row uses the reusable chromosome name `MT`:

```text
OG910.hifi.241127.v3mitohifi.emma102	MT	Circular-Chromosome	Mitochondrion
```

Packages are written beneath
`<outdir>/mitogenomes/<sample>/<assembly_prefix>/ena/package/`.
`<full_seqid>.package_metadata.json` is the single file the submission pipeline
reads: its `manifest` block holds the Webin genome-context keys under the names
Webin uses, and its `specimen` block holds the flatfile's source-feature facts.
`STUDY`, `SAMPLE` and `RUN_REF` are not in it, because the study, the sample and
the read submissions are registered downstream. [ena_package.md](ena_package.md)
documents the package and the hand-over to the submission pipeline.

The flatfile carries no `/locus_tag` on any feature. Locus tags are allocated and
injected by the downstream submission pipeline, so `table2asn` reports
`NO_LOCUS_TAGS` on every record and the pipeline treats that code as advisory.

`--ena_webin_validate` (on by default) format-validates every converted
flatfile with `ena-webin-cli -context sequence`. That is where the pipeline
stops, and passing it is what sets `submission_ready = true` on the validation
record. Genome-context validation needs a registered BioSample, so it belongs to
the downstream submission pipeline along with choosing and submitting a package.

ENA conversion runs automatically for samples with no table2asn `ERROR`/`REJECT` or discrepancy-report `FATAL`.
Warnings remain visible in MultiQC but do not block conversion. To add the credentialed production Webin check:

```bash
nextflow secrets set WEBIN_USERNAME 'Webin-XXXXXX'
nextflow secrets set WEBIN_PASSWORD

nextflow run main.nf \
  ... \
  --ena_webin_validate true
```

Every candidate in a run validates against one study, passed as `--ena_study`.
It has no default, so a run without it stops rather than validating against
whatever it fell back to.

The password is stored in Nextflow's secret store rather than a parameter or params file. The pipeline invokes
Webin with `-validate` only and never submits records. Individual Webin failures are written to the assembly's
`ena/validation/flatfile/` directory and do not terminate unrelated jobs. Only files under
`ena/validation/flatfile/validated/` have passed both the local table2asn gate and Webin validation.

### Standalone ENA conversion and validation

Use `ena.nf` when table2asn or EMBL outputs already exist and the assembly, annotation, and QC stages should not run.
The runner always performs Webin validation and therefore requires the two Webin secrets described above. Every row
validates against the run's single `--ena_study`, which is required.

To convert table2asn `.gbf` files and then validate them, create a CSV with these columns:

```csv
sample,mt_assembly_prefix,full_seqid,gbf,table2asn_status
OG1234,OG1234.mitohifi,OG1234.mitohifi.emma102,/path/OG1234.mitohifi.gbf,/path/OG1234.mitohifi.table2asn_status.tsv
```

`full_seqid` is optional and defaults to `mt_assembly_prefix`. It is what the validation record and
its database row are keyed on, so supply it whenever the flatfile was built from a named annotation;
without it the record's `annotation` column is empty.

The status file is mandatory. Only an exact table2asn `PASS` proceeds to conversion; failed or malformed statuses are
reported as `SKIP_TABLE2ASN` and do not stop other samples.

```bash
nextflow run ena.nf \
  -profile singularity \
  --ena_mode convert_validate \
  --ena_input ena_gbf_inputs.csv \
  --outdir results
```

Add `--sql_config /path/to/postgres.cfg` to either standalone command to upload the normalized
validation attempt to PostgreSQL after Webin finishes. Use `--skip_upload_results true` for an
explicitly file-only run. Nextflow never creates the database schema: apply the migrations first
(see [Database migrations](#database-migrations)).

`ena_validation_attempts` keeps one row per `(full_seqid, ena_study, validation_attempt)`.
Rerunning under the same attempt token overwrites that row rather than adding a new one, so retrying
a failed validation doesn't pile up history. Re-annotating an assembly is a different `full_seqid`
and therefore a new row, not an overwrite. Nothing freezes the row: submission state lives in the
separate submission pipeline, so the latest validation of a sequence is always the one on record.
Bump `--ena_validation_attempt` (see below) when you want a genuinely separate, independently
tracked attempt.

The PostgreSQL password remains in the protected SQL configuration file and is not written to ENA
manifests or normalized result records. Webin credentials continue to come only from Nextflow secrets.

To rerun Webin alone against existing compressed EMBL flat files, use:

```csv
sample,mt_assembly_prefix,full_seqid,embl
OG1234,OG1234.mitohifi,OG1234.mitohifi.emma102,/path/OG1234.mitohifi.embl.gz
```

```bash
nextflow run ena.nf \
  -profile singularity \
  --ena_mode validate \
  --ena_input ena_embl_inputs.csv \
  --outdir results
```

Relative input paths are resolved relative to the CSV location. Webin-only inputs first undergo gzip and EMBL
structure checks; only preflight-passing files reach Webin. Example CSVs are provided in
`assets/ena_gbf_samplesheet.csv` and `assets/ena_embl_samplesheet.csv`.

Webin failures are non-fatal and consequently cacheable. When using `-resume`, change the validation attempt token to
force Webin to check unchanged input again after a transient failure or credential change:

```bash
nextflow run ena.nf -resume \
  ... \
  --ena_validation_attempt 2026-07-15-retry1
```

The attempt token accepts letters, numbers, dots, underscores, and hyphens. It is recorded in each Webin status and
the combined `ena/ena_run_summary.tsv`. Standalone runs also create a MultiQC report containing the detailed gate
tables and the combined **ENA submission readiness** section.

## Mitogenome read depth

`MITOGENOME_COVERAGE` measures one depth number that means the same thing on every assembler and
platform: **mean per-base depth of the sample's own reads remapped to the assembly that goes to
annotation**. It writes `<assembly_prefix>.mito_depth.tsv` beside the assembly and populates
`mitogenome_data.mean_depth` in SQL.

This replaces three quantities that were never comparable with each other:

| Assembler | Legacy `avg_coverage` was | Now |
|-----------|---------------------------|-----|
| GetOrganelle | k-mer coverage off the assembly graph, roughly 0.2x true depth, over the reduced read set GetOrganelle selects by default | `mean_depth` |
| MitoHiFi | per-base depth of only the reads recruited by mapping to a related-species reference, so a divergent reference depressed it | `mean_depth` |
| Oatk | nothing at all (NULL) | `mean_depth` |

`avg_coverage` and `avg_base_coverage` are left untouched for provenance. Use `mean_depth` for any
cross-platform comparison, and `depth_method` to tell the two generations apart (`remap_full_v1`,
`not_measured`, or a `legacy_*` label applied by `sql/003_mitogenome_data_uniform_depth.sql`).

Two details worth knowing when reading the numbers:

- **Circular molecules are folded.** The assembly is doubled head to tail before mapping and the depth
  folded back, so origin-spanning reads are counted properly. Without this both ends of the linearised
  molecule show a false depth dip, which depresses the mean and inflates the CV.
- **NUMTs are filtered by gap-compressed identity**, so a read spanning a real control-region indel is
  kept while a diverged nuclear copy is rejected. MAPQ is deliberately ignored: on a doubled reference
  every read has two equally good placements, so MAPQ carries no information here.

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--skip_mitogenome_depth` | `false` | Skip the measurement; rows record `depth_method='not_measured'` |
| `--mitogenome_depth_min_identity_sr` | `0.95` | Minimum gap-compressed identity, Illumina/HiC |
| `--mitogenome_depth_min_identity_hifi` | `0.99` | Minimum gap-compressed identity, HiFi |
| `--mitogenome_depth_min_aligned_frac_sr` | `0.80` | Minimum aligned read fraction, Illumina/HiC (waived at reference ends) |
| `--mitogenome_depth_min_aligned_frac_hifi` | `0.70` | Minimum aligned read fraction, HiFi (a HiFi read can exceed the mitogenome length) |
| `--mitogenome_depth_subsample_fraction` | `0` | `0` uses every read. A value in (0,1) keeps a uniform random fraction and scales back up; saves mapping time but not decompression, so the full read set is the default |

Depth is measured only for the molecule that actually reaches annotation. Assemblies that failed, fell
below `--mitogenome_summary_min_length`, or are superseded GetOrganelle variants are recorded as
`not_measured` rather than remapped.

## Assembly summary QC thresholds

The pipeline writes `multiqc/mitogenome_assembly_summary_mqc.tsv` and includes it in MultiQC as
`Mitogenome assembly summary`. The parser uses these thresholds to populate `manual_review_reason`:

| Parameter | Default | Flag |
|-----------|---------|------|
| `--mitogenome_summary_min_mean_coverage` | `20` | `low_mean_coverage` (advisory on a complete assembly) |
| `--mitogenome_summary_max_coverage_cv` | `1.0` | `high_coverage_variability` (advisory on a complete assembly) |
| `--mitogenome_summary_min_length` | `10000` | `length_outside_expected_range` |
| `--mitogenome_summary_max_length` | `25000` | `length_outside_expected_range` |
| `--mitogenome_summary_expected_gene_count` | `37` | `missing_genes` when annotation-derived counts are lower |
| `--mitogenome_summary_expected_pcg_count` | `13` | `missing_protein_coding_genes` when the annotation-derived CDS count is lower (catches collapses that tRNAs mask in the total gene count) |

These thresholds are deliberately broad defaults for animal mitochondrial assemblies. Override them in
the command line or a params file when processing taxa with known compact, expanded, or unusual
mitogenomes.

The pipeline also writes one filtered MultiQC report per assembly, at
`mitogenomes/<sample>/<assembly_prefix>/<assembly_prefix>_multiqc_report.html`. Set
`--skip_per_sample_multiqc true` to generate only the cohort-level report.

## Samplesheet requirements (`--input`)

A valid CSV must match the schema in `assets/schema_input.json`:

- Required columns: `sample`, `fastq_1`, `sequencing_type`. `fastq_2` is optional.
- `sample` can be any non-empty identifier without whitespace.
- `sequencing_type` must be one of `ilmn`, `hifi`, or `hic`. This value controls which assembly
  subworkflow (GetOrganelle vs. MitoHiFi) is invoked and how metadata is derived.
- Optional metadata columns are accepted and propagated into `meta`, including
  `single_end`, `original_id`, `completion_date`, `date`, `assembly_prefix`,
  `nominal_species_id`, `invertebrates`, `class`, and `reference_species_id`.
- `class` (NCBI taxonomic class, e.g. `Actinopteri`, `Anthozoa`) and `invertebrates` together determine the
  per-sample mitochondrial genetic code used for annotation (e.g. Cnidaria → 4), the cox1 rotation panel, and
  which curated seed database a failed GetOrganelle first pass is reseeded from (see
  [`assets/refdb/README.md`](../assets/refdb/README.md)). In `--input_dir` mode these are
  resolved automatically from `nominal_species_id`: first from the OceanOmics `species` table, then from the NCBI
  taxdump (`--taxonkit_db_dir`) when the table has no match. In `--input` mode supply `class` yourself.
- **`class` is required in practice.** A missing or `unknown` class aborts the run before any assembly work, because
  it is indistinguishable from "vertebrate" downstream (genetic code 2, EMMA instead of MITOS2, the curated fish
  BLAST DB, the vertebrate tRNA completeness expectation). The error names every offending sample; fix them in the
  published `<outdir>/samplesheet/samplesheet.csv` (`taxonomy_resolution.tsv` alongside it shows where each lineage
  came from) and re-run with `--input`, or pass `--allow_unknown_taxonomy` to proceed with a warning.
- A blank `family`/`order` only warns: it stops the reference-divergence check grading beyond `NON_CONGENERIC`, so
  the `CROSS_ORDER` route to reference-free assembly cannot fire.
- When a sample has multiple libraries (e.g. several Illumina lanes), repeat the row with the same
  `sample` and `sequencing_type`. The pipeline concatenates the reads before downstream processing.

Example:

```csv title="oceanomics_samples.csv"
sample,fastq_1,fastq_2,sequencing_type,invertebrates,nominal_species_id
OG764,/data/OG764/OG764.ilmn.240716.R1.fastq.gz,/data/OG764/OG764.ilmn.240716.R2.fastq.gz,ilmn,false,Macruronus novaezelandiae
OG764,/data/OG764/OG764.ilmn.240717.R1.fastq.gz,/data/OG764/OG764.ilmn.240717.R2.fastq.gz,ilmn,false,Macruronus novaezelandiae
OG765_HICL,/data/OG765/OG765_HICL_S1_R1.fastq.gz,/data/OG765/OG765_HICL_S1_R2.fastq.gz,hic,false,Coelorinchus sp.
OG765,/pacbio/OG765_m84012_250101_s1.hifi_reads.fastq.gz,,hifi,true,Invertebrata sp.
```

When optional metadata fields such as `original_id` are present, they are propagated into `meta`
alongside `sample` so downstream modules can retain provenance. Do not insert whitespace in any
field.

Example sheets are available at `test_data/samplesheet.csv` (minimal) and `bin/samplesheet.csv`
(enriched real-world example).

## Directory mode (`--input_dir`)

Set `--input_dir "/path/to/OG*/**/*.fastq.gz"` to let the `CREATE_SAMPLESHEET_ENRICHED` module
generate the CSV for you. This mode requires `--sql_config`.

File names are parsed to determine sample identity, read pairing, sequencing mode, and collection
dates:

- Illumina: expects `<OGID>.ilmn.<YYMMDD>.R[12].fastq.gz` style names.
- PacBio HiFi: detects tokens such as `hifi_reads`, `hifi.reads`, or `.hifi` and records the date
  from the third underscore-delimited field.
- Hi-C: detects `HICL`/`HIC` in the file name and queries the sequencing date from the SQL database.

The generated CSV includes enriched metadata columns:
`sample,sequencing_type,single_end,original_id,completion_date,date,assembly_prefix,nominal_species_id,reference_species_id,class,family,order,invertebrates,genetic_code,fastq_1,fastq_2`.

It is published to `<outdir>/samplesheet/samplesheet.csv`, with `taxonomy_resolution.tsv` beside
it, and validated against `assets/schema_input.json`.

## External resource files

The pipeline orchestrates several helper modules that all require absolute file paths:

- **PostgreSQL configuration (`--sql_config`)** – INI file with credentials used by the enriched
  samplesheet builder (`--input_dir`) and SQL upload/QC modules. Example:
  ```ini
  [postgres]
  dbname = oceanomics
  user = analyst
  password = ******
  host = db.internal
  port = 5432
  ```
- **Curated BLAST database (`--curated_blast_db`)** – Provide the `*.fasta` (or BLAST database
  prefix) curated for OceanGenomes validations. Ensure `blastn` can resolve the files.
- **NCBI nt BLAST database (`--nt_blast_db`)** – Provide the `core_nt` (or equivalent) database
  prefix used for invertebrate samples. The pipeline selects this automatically when
  `invertebrates=true` in the samplesheet, and restricts the search to Metazoa
  (`-taxids 33208`) to keep a 1.09 Tbase database affordable.

  It must be a genuinely nt-scale database. Do **not** point this at the NCBI `mito`
  database: that is the RefSeq mitochondrion release (18,058 sequences), it holds almost
  none of the partial COI/16S submissions where invertebrate barcode evidence actually
  lives, and it produces empty BLAST output rather than an error. NOVA_260909_LA ran its
  whole invertebrate cohort that way — 4 of 23 samples reached the LCA, against 12 of 23
  on `core_nt`. Pass the database PREFIX (`.../db/core_nt`), not a file: the run fails
  fast if no `.nal` or `.nin` index is found beside it.
- **MITOS2 reference data (`--mitos_refdb`)** – RefSeq reference directory (parent of the version
  subdir named by `--mitos_refseq_ver`, default `refseq89m`) downloaded from Zenodo. Required to
  annotate invertebrate samples with MITOS2.
- **OatkDB profile-HMM (`--oatk_mito_db`)** – Only required when `--enable_oatk_fallback true`. The
  `<clade>_mito.fam` database plus its `.h3f/.h3i/.h3m/.h3p` index files; fetch with
  `bin/download_oatk_db.sh`.
- **Taxonomy caches (`--blast_db_dir`, `--taxonkit_db_dir`)** – Point these at a shared persistent
  filesystem. LCA prepares FishBase v24.07 and the NCBI taxdump once under
  `<taxonkit_db_dir>/lca_cache`, using a filesystem lock and atomic writes, and emits a checksum
  manifest before any per-gene LCA task can start. LCA tasks read that shared directory directly;
  the database is not copied or staged per task.
- **Submission template (`--template_sbt`)** – The `.sbt` template used by `table2asn` when
  packaging GenBank submission bundles.

Optional helpers include `--binddir` (Singularity bind mount root) and `--tempdir` (work directory
for large temporary files).

## Database migrations

Nextflow never creates or alters the database schema. The numbered migrations in `sql/` are applied
deliberately, before a run that uploads, with the migration runner:

```bash
python bin/apply_ena_migrations.py --config /path/to/postgres.cfg --check-only   # report pending
python bin/apply_ena_migrations.py --config /path/to/postgres.cfg                # apply pending
```

It applies the pending migrations from its list, in order, under an advisory lock, and audits the
schema before and after. Each one is recorded in `public.schema_migrations` (filename, sha256,
applied_at) and skipped once recorded, so a migration runs exactly once. `--baseline` records every
listed migration as applied *without running any of it*, for onboarding a database they were already
applied to; it refuses once the ledger is non-empty. Editing a migration after it has been applied
is an error (the database no longer matches the file), so the runner fails on checksum drift unless
you pass `--force`: write a new migration instead.

On Setonix there is no `psql` or `psycopg2` on the login nodes, so run it inside the psycopg2
container with `singularity exec <image> python bin/apply_ena_migrations.py ...`.

[`sql/README.md`](../sql/README.md) describes what each migration does and the ones that need care
(irreversible drops, table rebuilds, generated columns that must never be written).

## Skipping stages and reusing results

The main workflow supports coarse-grained skipping and reuse of pre-computed artefacts:

- `--skip_mitogenome_assembly_getorg` and `--skip_mitogenome_assembly_hifi` let you skip the
  assembly branches independently. Combine with:
  `--precomputed_mitogenome_assembly_fasta_getorg`,
  `--precomputed_mitogenome_assembly_log_getorg`,
  `--precomputed_mitogenome_assembly_fasta_hifi`,
  `--precomputed_mitogenome_assembly_log_hifi`.
- `--skip_getorganelle_reseed` – disable the automatic reseed pass so GetOrganelle only runs its
  first pass.
- `--skip_mitogenome_annotation` and `--blast_precomputed_annotations` – skip annotation and reuse
  existing results. These two flags behave differently from what their names suggest, so read
  [Skipping annotation](#skipping-annotation) before using them.
- `--skip_upload_results` – disable the SQL upload modules, and only those. Local QC is unaffected:
  annotation statistics, species validation, the QC gate, the QC summary and the held-samples report all
  still run, and the assembly summary still gets its annotation columns. ENA submission prep is unaffected
  too: it reads from the database without writing to it, so it runs whenever `--sql_config` is supplied.
  Use when working offline or on staging environments without database access.
- `--skip_ena_submission_prep` – disable table2asn, the ENA flat file and candidate packaging. Separate
  from `--skip_upload_results` because that flag no longer doubles as an off switch for this stage. Note
  that submission prep needs both `--sql_config` (it reads collection date, country and coordinates) and
  `--ena_study`; without a study it is skipped with a warning, since a study is what marks a run as doing
  submission work. `BUILD_SOURCE_MODIFIERS` fails rather than emitting an empty source-modifiers table for
  a sample whose `mitogenome_data` row does not exist yet.
- `--skip_hic_fastp` – skip fastp trimming of raw Hi-C reads (they are trimmed by default before
  GetOrganelle; Illumina reads bypass fastp because the draft-genome pipeline already trims them).

For precomputed assembly channels to map correctly back onto metadata, file basenames should encode
`<sample>.<sequencing_type>.<date>...` (for example:
`OG764.ilmn.240716.getorg1770.fasta`).

### Skipping annotation

Annotation, BLAST and LCA are three steps, but `--skip_mitogenome_annotation` does **not** skip
only the steps you have results for. It switches off all three together, and the steps after them
(QC, uploads, ENA prep) then read **all three** precomputed globs:

| Parameter | Default if you don't set it |
|---|---|
| `--precomputed_mitogenome_annotation_results` | `<outdir>/mitogenomes/*/*/annotation/*.{fa,fasta,gff,tbl,gb}` |
| `--precomputed_mitogenome_blast_results` | `<outdir>/mitogenomes/*/*/lca/blast.*filtered.tsv` |
| `--precomputed_mitogenome_lca_results` | `<outdir>/mitogenomes/*/*/lca/lca.*.tsv` |

Because every glob has a default, there is no such thing as "not providing" BLAST or LCA results.
If a glob matches no files, that step's results are simply **empty**. The pipeline does not error
and does not run BLAST instead. Species validation then sees no hits, and every sample is held as
`species_not_in_blast`.

To run BLAST + LCA on annotations you already have, add `--blast_precomputed_annotations true`:

| What you want | Flags | Annotation | BLAST + LCA |
|---|---|---|---|
| Normal run | (neither) | runs | runs |
| Reuse everything from an earlier run | `--skip_mitogenome_annotation true` | read from the annotation glob | read from the BLAST/LCA globs |
| Keep your annotations, redo BLAST + LCA (e.g. after hand curation) | `--skip_mitogenome_annotation true`<br>`--blast_precomputed_annotations true` | read from the annotation glob | **run** on those annotations; BLAST/LCA globs ignored |

In every case QC, uploads and ENA prep run afterwards as normal. `--blast_precomputed_annotations`
on its own, without `--skip_mitogenome_annotation`, does nothing (the run warns).

#### BLAST + LCA on existing annotations

Give the run its samplesheet as usual (`invertebrates` and `genetic_code` choose the BLAST database
and its arguments), the precomputed assemblies (the assembly rows and the QC gate need them), and the
annotation bundles:

```bash
SRC=/scratch/.../NOVA_260925_LA_mitogenomes/mitogenomes
nextflow run main.nf ... \
    --input samplesheet.csv \
    --skip_mitogenome_assembly_getorg true \
    --skip_mitogenome_assembly_hifi true \
    --precomputed_mitogenome_assembly_fasta_getorg "$SRC/*/*/mtdna/*.fasta" \
    --precomputed_mitogenome_assembly_log_getorg "$SRC/*/*/mtdna/*.get_org.log.txt" \
    --precomputed_mitogenome_circularity_evidence_getorg "$SRC/*/*/mtdna/*.getorg_check.tsv" \
    --skip_mitogenome_annotation true \
    --blast_precomputed_annotations true \
    --precomputed_mitogenome_annotation_results "$SRC/*/*/annotation/*.{fa,gff,tbl}"
```

- Glob only the bundle's top-level files, one annotation per assembly. Do not glob `cds/`: its
  `CO1.<prefix>.fa` names would be read as sample IDs, and the regions are re-extracted anyway. Two
  GFFs for one assembly (e.g. `mitos2110` and `mitos2110_mancurated`) is an error; narrow the glob.
- Each bundle's `cds/` and `proteins/` (which `ANNOTATION_STATS` needs for the PCG ORF check) and the
  CO1/12S/16S BLAST inputs are rebuilt from the GFF + FASTA, so an edited GFF is what gets QC'd and
  BLASTed. `annotation_rebuild/<annotation>.bundle_rebuild.tsv` records every protein-coding and rRNA
  gene and compares it with the bundle's own `cds/`: `match`, `trimmed` (EMMA's convention of dropping
  the stop codon; consistent with the GFF, so the bundle's `cds/` and `proteins/` are reused), `stale`
  (`cds/` disagrees with the GFF, which wins; the protein is re-translated), or `no_cds`.
- A region the GFF does not annotate is not BLASTed, exactly as on a normal run.

## Running with parameter files

Store frequently reused settings in a YAML/JSON file and pass it with `-params-file`:

```yaml title="params.yaml"
input: oceanomics_samples.csv
outdir: results
organelle_type: animal_mt
curated_blast_db: /databases/OceanGenomes.CuratedNT.fasta
nt_blast_db: /databases/blast/core_nt
sql_config: ~/.config/oceanomics/postgres.cfg
blast_db_dir: /scratch/shared/blast
taxonkit_db_dir: /scratch/shared/taxdump
template_sbt: ~/templates/oceanomics_submission.sbt
samplesheet_prefix: og_samples
translation_table: 2
skip_mitogenome_assembly_getorg: false
skip_mitogenome_assembly_hifi: false
```

```bash
nextflow run Minderoo-OceanOmics-Centre-UWA/Oceanomics-OceanGenomes-Mitogenomes -profile singularity -params-file params.yaml
```

Avoid using `-c` for parameter injection – reserve it for infrastructure or resource overrides.

## Execution profiles and `-resume`

All standard nf-core profiles are available (`docker`, `singularity`, `podman`, `conda`, `apptainer`,
`charliecloud`, `shifter`, `wave`). Pick the one that matches your environment; profiles can be
stacked (`-profile test,docker`).

Use `-resume` to restart a run from cached results. Nextflow will only reuse tasks whose inputs are
unchanged (including the SQL config, database paths, and FASTQ contents). Provide a run name
(`-resume <run-name>`) to target a specific execution recorded in `nextflow log`.

## Custom configuration

Resource requests follow the standard nf-core pattern (labels such as `process_low`, `process_high`
configured in `conf/base.config`). Customise them in a local config or via the nf-core/configs
repository if your institution requires bespoke defaults. To replace container images or inject
additional tool arguments, follow the guidance in the nf-core documentation:

- [Tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources)
- [Updating tool versions / containers](https://nf-co.re/docs/usage/configuration#updating-tool-versions)
- [Customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments)

## Stub runs and testing

Every process ships with a `stub` section so you can validate wiring and configuration quickly:

```bash
nextflow run Minderoo-OceanOmics-Centre-UWA/Oceanomics-OceanGenomes-Mitogenomes \
  -profile test,singularity \
  --input test_data/samplesheet.csv \
  --curated_blast_db test_data/blast_db \
  --nt_blast_db test_data/blast_db \
  --sql_config test_data/sql_config.txt \
  --template_sbt bin/template.sbt \
  --outdir stub_results \
  -stub-run
```

Stub runs generate placeholder files, which is useful for confirming connectivity (e.g. the SQL
config path) without launching long compute jobs.

## Output summary

See [output.md](output.md) for a stage-by-stage description of result files, including MultiQC
sections, SQL upload logs, and the GenBank and ENA package layout.
