# OceanOmics OceanGenomes Mitogenomes: Documentation

The OceanOmics OceanGenomes Mitogenomes documentation is split into the following pages:

- [Usage](usage.md)
  - How the pipeline works, how to run it (including `qc_only_from_annotations.nf` and `ena.nf`),
    the main command-line flags, and where the database schema lives. The
    full list, with defaults, is in `nextflow_schema.json` and `--help`.
- [Output](output.md)
  - The results directory, stage by stage, and how to interpret it.
- [ENA submission package](ena_package.md)
  - The package handed to the downstream ENA submission pipeline, how to build a Webin manifest
    from it, and the submission ledger.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)
