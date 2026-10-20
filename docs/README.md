# OceanOmics OceanGenomes Mitogenomes: Documentation

The OceanOmics OceanGenomes Mitogenomes documentation is split into the following pages:

- [Usage](usage.md)
  - How the pipeline works, how to run it (including `qc_only_from_annotations.nf` and `ena.nf`),
    every command-line flag, and applying database migrations.
- [Output](output.md)
  - The results directory, stage by stage, and how to interpret it.
- [ENA submission package](ena_package.md)
  - The package handed to the downstream ENA submission pipeline, how to build a Webin manifest
    from it, and the submission ledger.
- [Database migrations](../sql/README.md)
  - What each numbered migration in `sql/` does.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)
