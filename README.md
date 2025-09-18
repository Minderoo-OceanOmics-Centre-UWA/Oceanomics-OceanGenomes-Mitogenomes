# nf-core/oceangenomesmitogenomes

![nf-core/oceangenomesmitogenomes logo](docs/images/nf-core-oceangenomesmitogenomes_logo_light.png)

**nf-core/oceangenomesmitogenomes** is a modular [Nextflow DSL2](https://www.nextflow.io/) workflow that assembles,
annotates, validates, and packages mitochondrial genomes from PacBio HiFi, Illumina, and Hi-C sequencing runs. The
pipeline harmonises sample metadata, splits routing by sequencing mode, integrates with the OceanOmics PostgreSQL
database, and produces submission-ready GenBank artefacts alongside consolidated QC reporting.

## Introduction

OceanGenomes samples arrive from multiple sequencing modalities and require a consistent, auditable workflow to
convert raw reads into validated mitogenomes. This pipeline automates the complete journey:

- Samplesheet handling (`--input` CSV or `--input_dir` glob) with schema validation and automatic FASTQ discovery.
- Platform-aware assembly using GetOrganelle (short/Hi-C) and MitoHiFi (HiFi) with metadata-aware prefixes.
- Sanitisation, annotation (EMMA), BLAST validation, and lowest common ancestor (LCA) calling across key loci.
- Optional SQL upload modules that push assembly, annotation, and validation outputs into the OceanOmics database.
- QC gating and GenBank packaging (table2asn) for samples that satisfy species validation checks.
- MultiQC aggregation, software provenance capture, and full Nextflow run reports.

The workflow conforms to nf-core conventions, ships with reusable configuration profiles, and supports stub runs for
rapid wiring tests.

## Pipeline summary

- **Samplesheet preparation** – create or validate CSV input, derive sequencing dates, and enrich metadata for routing.
- **Mitogenome assembly** – concatenate lane-level FASTQs and run GetOrganelle or MitoHiFi depending on sequencing type.
- **Sanitise and annotate** – normalise FASTA headers, run EMMA, and export gene- and protein-level FASTA collections.
- **BLAST + LCA** – download taxonomy caches on-demand, filter BLAST hits for CO1/12S/16S, and compute per-sample LCA.
- **Database integration** – push assembly metrics, annotation stats, and BLAST/LCA summaries to PostgreSQL (optional).
- **QC and packaging** – evaluate validation criteria, build GenBank-ready bundles, run `table2asn`, and collect diagnostics.
- **Reporting** – collate MultiQC inputs, generate run-level reports, and write `pipeline_info/` provenance files.

See `docs/images` for workflow logos and refer to `docs/output.md` for a stage-by-stage file manifest.

## Quick start

1. Install Nextflow >= 24.10.0 and a container backend (Singularity, Docker, Podman, or Conda/Mamba).
2. Pull the pipeline or clone this repository: `nextflow pull nf-core/oceangenomesmitogenomes`.
3. Launch a run with an execution profile and an input definition.

### Samplesheet mode (`--input`)

```bash
nextflow run nf-core/oceangenomesmitogenomes \
  -profile singularity \            # or docker / podman / conda / mamba
  --input path/to/oceanomics_samples.csv \
  --outdir results \
  --organelle_type animal_mt \
  --curated_blast_db /path/to/OceanGenomes.CuratedNT.fasta \
  --sql_config ~/postgresql_details/oceanomics.cfg \
  --blast_db_dir /scratch/shared/blast \
  --taxonkit_db_dir /scratch/shared/taxdump \
  --template_sbt ~/templates/oceanomics_submission.sbt
```

### Directory mode (`--input_dir`)

```bash
nextflow run nf-core/oceangenomesmitogenomes \
  -profile singularity \
  --input_dir "/data/OG*/**/*.fastq.gz" \
  --samplesheet_prefix og_run \
  --outdir results \
  --organelle_type animal_mt \
  --curated_blast_db /path/to/OceanGenomes.CuratedNT.fasta \
  --sql_config ~/postgresql_details/oceanomics.cfg \
  --blast_db_dir /scratch/shared/blast \
  --taxonkit_db_dir /scratch/shared/taxdump \
  --template_sbt ~/templates/oceanomics_submission.sbt
```

### Stub run for wiring/tests

```bash
nextflow run nf-core/oceangenomesmitogenomes \
  -profile test,singularity \  # uses bundled test data
  --input assets/samplesheet.csv \
  --curated_blast_db test_data/blast_db \
  --sql_config test_data/sql_config.txt \
  --template_sbt test_data/template.sbt \
  --outdir stub_results \
  -stub-run
```

### Pawsey example

A reference launch script tailored for Pawsey (Setonix) is provided in `nextflow_run.sh`. It loads the recommended
modules, applies `pawsey_profile.config`, stages cached databases under `/scratch`, and toggles skip flags. Use it as a
starting point for production-scale runs on that system.

## Documentation

Comprehensive usage and output documentation is maintained under `docs/`:

- `docs/usage.md` – parameter descriptions, samplesheet schema, external resource requirements, and execution advice.
- `docs/output.md` – per-stage file inventories, MultiQC sections, and guidance on interpreting database/QC artefacts.
- `docs/README.md` – entry point with links to the above.

For general nf-core best practices (custom configs, resource tuning, module overrides), see <https://nf-co.re/docs>.

## Key parameters and resources

- `--input` / `--input_dir` – mutually exclusive input modes; both validate against `assets/schema_input.json`.
- `--organelle_type` – passed to GetOrganelle (e.g. `animal_mt`, `embryophyta_mt`).
- `--curated_blast_db` – absolute path to the OceanGenomes-curated BLAST database.
- `--sql_config` – INI credentials for PostgreSQL-backed modules (`SPECIES_QUERY`, upload steps, QC gating).
- `--blast_db_dir`, `--taxonkit_db_dir` – scratch/persistent directories for caching taxonomy resources.
- `--template_sbt` – GenBank submission template consumed by `table2asn`.
- `--skip_*` and `--precomputed_*` – coarse-grained switches for reusing existing assemblies, annotations, or validation
  outputs. See `nextflow_schema.json` for the full inventory.

All available parameters (including hidden/advanced options) are documented in the autogenerated help:

```bash
nextflow run nf-core/oceangenomesmitogenomes --help
```

## Configuration profiles

This repository ships with standard nf-core profiles (`test`, `test_full`, `docker`, `singularity`, `podman`,
`conda`, `apptainer`, `wave`, etc.) defined in `nextflow.config`. Additional helpers include:

- `-profile pawsey` (via `pawsey_profile.config`) – Slurm + Singularity defaults for Setonix.
- `-stub-run` – execute stub sections for rapid validation of channel wiring and configuration.
- `-profile debug` – enable verbose logging, hash dumps, and DAG validation.

Profiles can be combined (e.g. `-profile test,singularity`). Use `-params-file` to store recurring parameters in YAML or
JSON and keep `-c` for infrastructure overrides.

## Results and reporting

Final results are published under `--outdir` with predictable subdirectories:

- `mitogenomes/` – assemblies, annotations, BLAST/LCA outputs, and GenBank packaging artefacts per sample.
- `species_validation/` and `sql_uploaded_data/` – per-sample SQL upload logs and QC summaries.
- `multiqc/` – consolidated HTML report plus supporting data.
- `pipeline_info/` – Nextflow reports, parameter snapshots, software versions, and validated samplesheet copies.

Refer to `docs/output.md` for exhaustive listings.

## Testing

- `-profile test,<container>` exercises the pipeline with minimal bundled data.
- `tests/` and `nf-test.config` provide scaffolding for [`nf-test`](https://github.com/nf-core/nf-test) based unit tests.
- Continuous integration via GitHub Actions can be adapted from the nf-core template to automate linting and `nf-test`.

## Credits

Developed by Tyler Peirce and the OceanGenomes/OceanOmics bioinformatics team. See `CHANGELOG.md` for release history and
`modules/` for third-party module acknowledgements. Logo assets are stored under `docs/images/` and `assets/`.

## Contributions

Contributions are welcome! Please review `CODE_OF_CONDUCT.md`, open an issue for major changes, and follow the nf-core
module and pipeline development guidelines when submitting pull requests. Linting, documentation updates, and tests are
highly appreciated.

## Citations

If you use nf-core/oceangenomesmitogenomes in your work, please cite:

1. The nf-core framework: `Ewels et al., Nature Biotechnology 38, 276-278 (2020)` – <https://doi.org/10.1038/s41587-020-0439-x>
2. Nextflow: `Di Tommaso et al., Nature Biotechnology 35, 316-319 (2017)` – <https://doi.org/10.1038/nbt.3820>
3. The specific tools referenced in `CITATIONS.md` (EMMA, GetOrganelle, MitoHiFi, BLAST, TaxonKit, MultiQC, etc.).

A machine-readable bibliography is provided in `CITATIONS.md` for convenience.
