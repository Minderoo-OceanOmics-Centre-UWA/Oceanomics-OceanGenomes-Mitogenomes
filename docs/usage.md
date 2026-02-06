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

Nextflow ≥ 24.10.x is required. Always combine the pipeline with an execution profile (`-profile`)
that matches your container/conda environment.

## Essential parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--input` | ✔ (or `--input_dir`) | Path to a CSV samplesheet that matches `assets/schema_input.json`. |
| `--input_dir` | ✔ (or `--input`) | Glob pointing at FASTQ files when you want the pipeline to build the samplesheet automatically. |
| `--outdir` | ✔ | Destination for all published results. |
| `--organelle_type` | ✔ | Organelle label passed to GetOrganelle (e.g. `animal_mt`). |
| `--curated_blast_db` | ✔ | NCBI-format BLAST database used for species validation. Provide an absolute path. |
| `--nt_blast_db` | Conditional | NCBI nt BLAST database used when samples are marked `invertebrates=true`. |
| `--sql_config` | ✔* | INI file with `[postgres] dbname,user,password,host,port`. Required whenever Hi-C samples or species validation/database upload steps run. |
| `--blast_db_dir` | ✔ | Directory used to cache the downloaded `taxdb.*` files; re-use between runs to avoid repeated downloads. |
| `--taxonkit_db_dir` | ✔ | Directory used to cache the NCBI taxdump for TaxonKit. |
| `--template_sbt` | ✔ | Submission template passed to `table2asn` when packaging GenBank artefacts. |
| `--samplesheet_prefix` | Optional | Prefix for the auto-generated samplesheet when using `--input_dir`. Defaults to `samplesheet`. |
| `--translation_table` | Optional | Override the mitochondrial translation table (defaults to vertebrate code `2`). |

`--sql_config` can be omitted only when **all** modules that touch the OceanOmics PostgreSQL
instance are skipped (for example, by running `-stub-run` or disabling validation/upload stages).
Otherwise, the pipeline will raise an error early on.

## Samplesheet requirements (`--input`)

A valid CSV must match the schema in `assets/schema_input.json`:

- Required columns: `sample`, `fastq_1`, `sequencing_type`. `fastq_2` is optional.
- `sample` must follow the OceanGenomes convention (`OG` followed by digits, optional suffixes such
  as `OG820M-1`). The schema enforces this naming pattern.
- `sequencing_type` must be one of `ilmn`, `hifi`, or `hic`. This value controls which assembly
  subworkflow (GetOrganelle vs. MitoHiFi) is invoked and how metadata is derived.
- When a sample has multiple libraries (e.g. several Illumina lanes), repeat the row with the same
  `sample` and `sequencing_type`. The pipeline concatenates the reads before downstream processing.

Example:

```csv title="oceanomics_samples.csv"
sample,fastq_1,fastq_2,sequencing_type
OG764.ilmn.240716,/data/OG764/OG764.ilmn.240716.R1.fastq.gz,/data/OG764/OG764.ilmn.240716.R2.fastq.gz,ilmn
OG764.ilmn.240716,/data/OG764/OG764.ilmn.240717.R1.fastq.gz,/data/OG764/OG764.ilmn.240717.R2.fastq.gz,ilmn
OG765.hic,/data/OG765/OG765_HICL_S1_R1.fastq.gz,/data/OG765/OG765_HICL_S1_R2.fastq.gz,hic
OG765.hifi,/pacbio/OG765_m84012_250101_s1.hifi_reads.fastq.gz,,hifi
```

To keep the schema strict while preserving the original identifiers, the pipeline stores the
provided ID alongside a cleaned `meta.id`. Do not insert whitespace in any field.

An example sheet lives at `assets/samplesheet.csv`.

## Directory mode (`--input_dir`)

Set `--input_dir "/path/to/OG*/**/*.fastq.gz"` to let the `CREATE_SAMPLESHEET` module generate the
CSV for you. File names are parsed to determine sample identity, read pairing, sequencing mode, and
collection dates:

- Illumina: expects `<OGID>.ilmn.<YYMMDD>.R[12].fastq.gz` style names.
- PacBio HiFi: detects tokens such as `hifi_reads`, `hifi.reads`, or `.hifi` and records the date
  from the third underscore-delimited field.
- Hi-C: detects `HICL`/`HIC` in the file name and queries the sequencing date from the SQL database
  via `HIC_DATE_QUERY`.

Use `--samplesheet_prefix` to control the generated filename (defaults to `samplesheet.csv`). The
file is written to the working directory and validated against the same schema as manual inputs.

## External resource files

The pipeline orchestrates several helper modules that all require absolute file paths:

- **PostgreSQL configuration (`--sql_config`)** – INI file with credentials used by
  `SPECIES_QUERY`, `HIC_DATE_QUERY`, and the upload modules. Example:
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

Optional helpers include `--binddir` (Singularity bind mount root), `--tempdir` (work directory for
large temporary files), and `--refresh-modules` (force module cache refresh when running on the
OceanOmics infrastructure scripts).

## Skipping stages and reusing results

The main workflow supports coarse-grained skipping and reuse of pre-computed artefacts:

- `--skip_mitogenome_assembly` – bypass both GetOrganelle and MitoHiFi. Combine with
  `--precomputed_mitogenome_assembly_fasta` and
  `--precomputed_mitogenome_assembly_log` to supply existing assemblies/logs.
- `--skip_mitogenome_annotation` – skip EMMA / BLAST / LCA; optional inputs
  `--precomputed_mitogenome_annotation_results`, `--precomputed_mitogenome_blast_results`, and
  `--precomputed_mitogenome_lca_results` allow you to feed downstream modules.
- `--skip_upload_results` – disable SQL upload modules and QC gating. Use when working offline or on
  staging environments without database access.

Refer to `nextflow_schema.json` for the complete list of `skip_*` and `precomputed_*` parameters.

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
  --input assets/samplesheet.csv \
  --curated_blast_db test_data/blast_db \
  --nt_blast_db test_data/blast_db \
  --sql_config test_data/sql_config.txt \
  --template_sbt test_data/template.sbt \
  --outdir stub_results \
  -stub-run
```

Stub runs generate placeholder files, which is useful for confirming connectivity (e.g. the SQL
config path) without launching long compute jobs.

## Output summary

See `docs/output.md` for a detailed, stage-by-stage description of result files, including MultiQC
sections, SQL upload logs, and the GenBank submission bundle layout.
