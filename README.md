## oceangenomes-mitogenomes

End-to-end workflow for mitochondrial genome assembly, annotation, species validation and submission-ready QC. Supports PacBio HiFi, Illumina, and Hi-C datasets with automatic routing and consistent sample metadata handling.

Workflow overview
- Input handling: use a CSV samplesheet or point to a FASTQ glob; a samplesheet will be generated if needed.
- Sequencing routing: branch by `sequencing_type` (hifi → MitoHiFi; ilmn/hic → GetOrganelle). Compute a consistent `assembly_prefix` per sample.
- Annotation and validation: annotate mitogenomes (EMMA), run BLAST against a curated database, compute LCA for species validation.
- Optional database steps: query PostgreSQL for Hi-C sequencing date and nominal species; optionally upload LCA/BLAST results.
- QC and submission: extract genes, translate proteins, format files and generate GenBank submission artifacts (e.g., `.sqn`, `.val`, `.gbf`).
- Reporting: aggregate MultiQC and collate software versions.

Inputs
- Samplesheet mode: `--input` CSV matching `assets/schema_input.json` with columns: `sample,fastq_1,fastq_2,sequencing_type`.
- Directory mode: `--input_dir` glob for FASTQs and `--samplesheet_prefix` to auto-generate a CSV.

Key parameters
- `--input`: Path to samplesheet CSV.
- `--input_dir`: Glob for FASTQ inputs; triggers samplesheet generation.
- `--samplesheet_prefix`: Filename prefix when generating the samplesheet.
- `--outdir`: Output directory root.
- `--organelle_type`: Organelle target for assembly (e.g., `animal_mt`).
- `--sql_config` (optional): INI with `[postgres] dbname,user,password,host,port` for DB steps.
- `--curated_blast_db` (optional): Path to curated BLAST database for annotation/validation.
- Skips/precomputed: `skip_*` and `precomputed_*` toggles available in `nextflow.config` to reuse results.

Quick start
- Nextflow ≥ 24.10.5 (DSL2)
- Profiles: docker / singularity / conda

Samplesheet input
```
nextflow run . \
  -profile <docker|singularity|conda> \
  --input /path/to/samples.csv \
  --outdir results \
  --organelle_type animal_mt
```

Directory input (auto samplesheet)
```
nextflow run . \
  -profile <docker|singularity|conda> \
  --input_dir "/data/*.fastq.gz" \
  --samplesheet_prefix samples \
  --outdir results \
  --organelle_type animal_mt
```

Stub run
- All modules implement `stub:` blocks for fast dry-runs:
```
nextflow run . -stub-run -profile test,<docker|singularity|conda> \
  --input /path/to/samples.csv \
  --outdir stub_results
```

Outputs (high-level)
- Assembly: FASTA plus logs per sample (GetOrganelle and/or MitoHiFi).
- Annotation and validation: EMMA outputs, BLAST filtered results, LCA summaries, species decision.
- Submission QC: directories with translated proteins and generated submission files (`*.sqn`, `*.val`, `*.gbf`).
- Reports: MultiQC HTML and a collated `pipeline_info` directory with versions and traces.

Notes
- See `nextflow_schema.json` for validated parameters and defaults; profiles are in `nextflow.config`.
- DB-related modules no-op under `-stub-run` and emit placeholder outputs.

Credits
- Developed by Tyler Peirce and collaborators.
