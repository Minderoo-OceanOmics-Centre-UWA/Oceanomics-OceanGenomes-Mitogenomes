# nf-core/oceangenomesmitogenomes: Usage

The nf-core/oceangenomesmitogenomes workflow assembles, annotates, validates, and packages
mitochondrial genomes from PacBio HiFi, Illumina, and Hi-C sequencing runs. This page describes
pipeline-specific input requirements and the parameters you must provide to execute a full run on
OceanOmics data holdings.

> For general guidance on Nextflow execution, configuration profiles, and infrastructure tuning,
> refer to the nf-core documentation: <https://nf-co.re/docs/usage>

## Quick start

```bash
nextflow run nf-core/oceangenomesmitogenomes \
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
| `--translation_table` | Optional | Fallback mitochondrial genetic code for gene translation when no per-sample code is available (defaults to `2`). In QC-only mode the taxonomic class is not re-derived, so this value is used directly. |

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
| `--nt_blast_db` | Conditional | NCBI nt BLAST database used when samples are marked `invertebrates=true`. |
| `--sql_config` | Conditional | INI file with `[postgres] dbname,user,password,host,port`. Required for `--input_dir`; optional for `--input` if upload/QC is skipped. |
| `--blast_db_dir` | ✔ | Directory used to cache the downloaded `taxdb.*` files; re-use between runs to avoid repeated downloads. |
| `--taxonkit_db_dir` | ✔ | Directory used to cache the NCBI taxdump for TaxonKit. |
| `--template_sbt` | ✔ | Submission template passed to `table2asn` when packaging GenBank artefacts. |
| `--ena_webin_validate` | Optional | Run production Webin-CLI validation for ENA flat files that pass table2asn and conversion checks (default `false`). |
| `--ena_study` | With Webin validation | ENA study accession or alias written to each targeted-sequence manifest. |
| `--samplesheet_prefix` | Optional | Reserved for generated samplesheet naming in wrapper scripts. |
| `--getorganelle_genedb_min_genes` | Optional | Minimum genes a reference must yield to build the reseed custom gene database (default `10`). Below this, the sample keeps its first-pass GetOrganelle assembly instead of reseeding. |
| `--translation_table` | Optional | Mitochondrial genetic code for vertebrate/unresolved samples (default `2`). Invertebrate codes are derived per-sample from the taxonomic `class` column (Cnidaria → 4, echinoderms/flatworms → 9, other invertebrates → 4), so this no longer forces a single code across the whole run. |

`--input_dir` mode requires `--sql_config` because the enriched samplesheet generator queries
OceanOmics metadata. When running with `--input`, you can omit `--sql_config`, but upload/QC stages
are then skipped with a warning.

### ENA Webin validation

ENA conversion runs automatically for samples with no table2asn `ERROR`/`REJECT` or discrepancy-report `FATAL`.
Warnings remain visible in MultiQC but do not block conversion. To add the credentialed production Webin check:

```bash
nextflow secrets set WEBIN_USERNAME 'Webin-XXXXXX'
nextflow secrets set WEBIN_PASSWORD

nextflow run main.nf \
  ... \
  --ena_webin_validate true \
  --ena_study PRJEB12345
```

The password is stored in Nextflow's secret store rather than a parameter or params file. The pipeline invokes
Webin with `-validate` only and never submits records. Individual Webin failures are written to the sample's
`genbank/ena/` directory and do not terminate unrelated jobs. Only files under `ena/validated/` have passed both
the local table2asn gate and Webin validation.

### Standalone ENA conversion and validation

Use `ena.nf` when table2asn or EMBL outputs already exist and the assembly, annotation, and QC stages should not run.
The runner always performs Webin validation and therefore requires `--ena_study` and the two Webin secrets described
above.

To convert table2asn `.gbf` files and then validate them, create a CSV with these columns:

```csv
sample,mt_assembly_prefix,gbf,table2asn_status
OG1234,OG1234.mitohifi,/path/OG1234.mitohifi.gbf,/path/OG1234.mitohifi.table2asn_status.tsv
```

The status file is mandatory. Only an exact table2asn `PASS` proceeds to conversion; failed or malformed statuses are
reported as `SKIP_TABLE2ASN` and do not stop other samples.

```bash
nextflow run ena.nf \
  -profile singularity \
  --ena_mode convert_validate \
  --ena_input ena_gbf_inputs.csv \
  --ena_study PRJEB123456 \
  --outdir results
```

Add `--sql_config /path/to/postgres.cfg` to either standalone command to append the normalized
validation attempt to PostgreSQL after Webin finishes. Use `--skip_upload_results true` for an
explicitly file-only run. The database schema is never created by Nextflow; apply the migration
deliberately before enabling this upload, for example:

```bash
psql --dbname oceanomics --file sql/001_create_ena_validation_attempts.sql
```

The PostgreSQL password remains in the protected SQL configuration file and is not written to ENA
manifests or normalized result records. Webin credentials continue to come only from Nextflow secrets.

To rerun Webin alone against existing compressed EMBL flat files, use:

```csv
sample,mt_assembly_prefix,embl
OG1234,OG1234.mitohifi,/path/OG1234.mitohifi.embl.gz
```

```bash
nextflow run ena.nf \
  -profile singularity \
  --ena_mode validate \
  --ena_input ena_embl_inputs.csv \
  --ena_study PRJEB123456 \
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

## Assembly summary QC thresholds

The pipeline writes `multiqc/mitogenome_assembly_summary_mqc.tsv` and includes it in MultiQC as
`Mitogenome assembly summary`. The parser uses these thresholds to populate `manual_review_reason`:

| Parameter | Default | Flag |
|-----------|---------|------|
| `--mitogenome_summary_min_mean_coverage` | `20` | `low_mean_coverage` |
| `--mitogenome_summary_max_coverage_cv` | `1.0` | `high_coverage_variability` |
| `--mitogenome_summary_min_length` | `10000` | `length_outside_expected_range` |
| `--mitogenome_summary_max_length` | `25000` | `length_outside_expected_range` |
| `--mitogenome_summary_expected_gene_count` | `37` | `missing_genes` when annotation-derived counts are lower |

These thresholds are deliberately broad defaults for animal mitochondrial assemblies. Override them in
the command line or a params file when processing taxa with known compact, expanded, or unusual
mitogenomes.

The pipeline also writes one filtered MultiQC report per detected sample under
`multiqc/per_sample/<sample>/<sample>_multiqc_report.html`. Set `--skip_per_sample_multiqc true`
to generate only the cohort-level report.

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
  per-sample mitochondrial genetic code used for annotation (e.g. Cnidaria → 4). In `--input_dir` mode these are
  resolved automatically from `nominal_species_id`; in `--input` mode supply `class` for correct invertebrate codes.
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
`sample,sequencing_type,single_end,original_id,completion_date,date,assembly_prefix,nominal_species_id,invertebrates,fastq_1,fastq_2`.

The current implementation writes `samplesheet.csv` and validates it against
`assets/schema_input.json`.

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
  `invertebrates=true` in the samplesheet.
- **Taxonomy caches (`--blast_db_dir`, `--taxonkit_db_dir`)** – Large downloads (`taxdb.*`, NCBI
  `taxdump.tar.gz`) are stored under these directories using `storeDir`. Point them to a shared or
  persistent filesystem to avoid repeated downloads.
- **Submission template (`--template_sbt`)** – The `.sbt` template used by `table2asn` when
  packaging GenBank submission bundles.

Optional helpers include `--binddir` (Singularity bind mount root) and `--tempdir` (work directory
for large temporary files).

## Skipping stages and reusing results

The main workflow supports coarse-grained skipping and reuse of pre-computed artefacts:

- `--skip_mitogenome_assembly_getorg` and `--skip_mitogenome_assembly_hifi` let you skip the
  assembly branches independently. Combine with:
  `--precomputed_mitogenome_assembly_fasta_getorg`,
  `--precomputed_mitogenome_assembly_log_getorg`,
  `--precomputed_mitogenome_assembly_fasta_hifi`,
  `--precomputed_mitogenome_assembly_log_hifi`.
- `--skip_mitogenome_annotation` – skip EMMA / BLAST / LCA; optional inputs
  `--precomputed_mitogenome_annotation_results`, `--precomputed_mitogenome_blast_results`, and
  `--precomputed_mitogenome_lca_results` allow you to feed downstream modules.
- `--skip_upload_results` – disable SQL upload modules and QC gating. Use when working offline or on
  staging environments without database access.

For precomputed assembly channels to map correctly back onto metadata, file basenames should encode
`<sample>.<sequencing_type>.<date>...` (for example:
`OG764.ilmn.240716.getorg1770.fasta`).

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
nextflow run nf-core/oceangenomesmitogenomes -profile singularity -params-file params.yaml
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
nextflow run nf-core/oceangenomesmitogenomes \
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

See `docs/output.md` for a detailed, stage-by-stage description of result files, including MultiQC
sections, SQL upload logs, and the GenBank submission bundle layout.
