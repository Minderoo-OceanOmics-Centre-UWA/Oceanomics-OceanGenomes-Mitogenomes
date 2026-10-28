# OceanOmics OceanGenomes Mitogenomes: ENA submission package

This pipeline builds and format-validates an ENA submission package for every viable assembly and
annotation, but it does not submit anything. Submission is done by a separate pipeline,
[ENA-mito-genomes](https://github.com/Minderoo-OceanOmics-Centre-UWA/ENA-mito-genomes). This page
describes the package that is handed over, how to build a Webin manifest from it, and the database
tables the two pipelines share.

## Responsibilities

Specimens sit under embargo for an arbitrary period after their mitogenome is assembled, so
submission cannot be a stage of the assembly pipeline: the two events are months apart. The hand-over
is a directory of files.

| | This pipeline | ENA-mito-genomes |
|---|---|---|
| Builds the EMBL flat file and chromosome list | yes | no |
| Format-validates the flat file (`webin-cli -context sequence -validate`) | yes | no |
| Builds the submission manifest | no | yes |
| Allocates and injects locus tags | no | yes |
| Decides which assembly represents a specimen | no | yes |
| Applies the embargo gate | no | yes |
| Calls `webin-cli -context genome -submit` | no | yes |
| Records the resulting accessions | defines the table | writes it |

This pipeline only ever runs `webin-cli -context sequence -validate`, against the flat file alone.
Genome-context validation cannot run here at all: it resolves `SAMPLE` against ENA's submission
sample service before it checks anything else, and the sample is registered downstream.

## Where packages are

Each run writes packages to `<outdir>/mitogenomes/<og_id>/<assembly_prefix>/ena/package/`. The
backup script (`backup_scripts/backup.sh`) copies them to:

```text
s3://ocom.oceangenomes/analysed-data/mitogenomes/curated/<og_id>/<assembly_prefix>/
```

Each package describes itself; reading one needs no database join.

## Package contents

Every file is named on `full_seqid` (the assembly prefix plus the annotation version, e.g.
`OG51.ilmn.240313.getorg1770.emma102`), so the whole package shares one stem.

| File | Purpose |
|---|---|
| `<full_seqid>.embl.gz` | the flat file; Webin `FLATFILE` input |
| `<full_seqid>.chromosome_list.tsv.gz` | Webin `CHROMOSOME_LIST` input |
| `<full_seqid>.package_metadata.json` | everything this pipeline knows about the submission |
| `<full_seqid>.tbl` | NCBI feature table, kept for provenance |
| `<full_seqid>.fa` | sequence, for collaborators |
| `<full_seqid>.gff` | annotation as the annotator produced it, for collaborators |
| `<full_seqid>.genes.fa` | per-gene sequences (CDS, tRNA, rRNA), for collaborators |
| `checksums.sha256` | sha256 of every other file in the directory |

The `.fa`, `.gff` and `.genes.fa` are **not** submission inputs.

The chromosome list is a single row:

```text
<full_seqid>	MT	Circular-Chromosome	Mitochondrion
```

## Building the manifest

The package contains no manifest file, because three of the Webin genome-context keys are values
only the submission pipeline holds:

- **`STUDY`**: the BioProject the assembly is submitted into.
- **`SAMPLE`**: the BioSample registered for the specimen.
- **`RUN_REF`**: the raw-read submissions.

The `manifest` block of `package_metadata.json` holds the other nine keys, already under the names
Webin uses. Building the manifest is copying that block out as tab-separated lines and adding those
three. Injecting locus tags is the only change the package itself needs before `-submit`.

The `manifest` block for OG51:

```json
"manifest": {
  "ASSEMBLYNAME": "OG51.ilmn.240313.getorg1770.emma102",
  "ASSEMBLY_TYPE": "clone or isolate",
  "CHROMOSOME_LIST": "OG51.ilmn.240313.getorg1770.emma102.chromosome_list.tsv.gz",
  "COVERAGE": "1698.16",
  "DESCRIPTION": "Ophthalmolepis lineolata mitochondrial genome",
  "FLATFILE": "OG51.ilmn.240313.getorg1770.emma102.embl.gz",
  "MOLECULETYPE": "genomic DNA",
  "PLATFORM": "ILLUMINA",
  "PROGRAM": "GetOrganelle 1.7.7.0"
}
```

`FLATFILE` and `CHROMOSOME_LIST` are bare file names, relative to the package directory.

A key whose value would not validate is left out rather than written wrong, so treat the block as
authoritative about what it contains, not as a fixed set of nine. A missing key means the pipeline
had nothing trustworthy to put there.

| Manifest key | Source |
|---|---|
| `STUDY`, `SAMPLE`, `RUN_REF` | the submission pipeline |
| `ASSEMBLYNAME` | `full_seqid` |
| `ASSEMBLY_TYPE` | fixed: `clone or isolate` |
| `COVERAGE` | `mean_depth`, rendered with `%g` |
| `PROGRAM`, `PLATFORM` | the assembler and sequencing platform, in Webin vocabulary |
| `MOLECULETYPE` | fixed: `genomic DNA` |
| `DESCRIPTION` | `<scientific_name> mitochondrial genome` |
| `FLATFILE`, `CHROMOSOME_LIST` | the package file names |

The manifest and the flat file are separate webin-cli inputs, checked separately. What the flat file
carries does not license leaving anything out of the manifest.

### Flat file and manifest

- **`full_seqid` is the `ASSEMBLYNAME`**, the file name stem and the flat file `ID` line, so all three
  agree by construction. It includes the annotation version, so a re-annotation is a new assembly
  rather than a collision with one already submitted.
- **The `specimen` block** carries the qualifiers of the flat file's `source` feature, read back out
  of the packaged flat file, so the JSON and the flat file cannot disagree. A qualifier with no valid
  value is **absent, not empty** (OG51 has no `/lat_lon`, so there is no `lat_lon` key), and
  `organism` is normalised (`Genus sp` reads `Genus sp.`).

  ```json
  "specimen": {
    "collection_date": "29-Mar-2023",
    "geo_loc_name": "Australia: WA, New Year Island",
    "isolate": "OG51",
    "mol_type": "genomic DNA",
    "organelle": "mitochondrion",
    "organism": "Ophthalmolepis lineolata",
    "tissue_type": "Gills"
  }
  ```

- **The flat file carries no study, sample, run or coverage.** What it does carry:

  ```text
  ID   OG51.ilmn.240313.getorg1770.emma102; SV 1; circular; genomic DNA; STD; UNC; 16631 BP.
  XX
  AC   ;
  XX
  AC * _OG51.ilmn.240313.getorg1770.emma102
  XX
  DE   Ophthalmolepis lineolata mitochondrion, complete genome
  XX
  OS   Ophthalmolepis lineolata
  XX
  CC   ##Assembly-Data-START##
  CC   Assembly Method       :: GetOrganelle v.1.7.7.0
  CC   Sequencing Technology :: Illumina
  CC   ##Assembly-Data-END##
  ...
  FT   source          1..16631
  FT                   /organism="Ophthalmolepis lineolata"
  FT                   /organelle="mitochondrion"
  FT                   /mol_type="genomic DNA"
  FT                   /isolate="OG51"
  FT                   /tissue_type="Gills"
  FT                   /geo_loc_name="Australia: WA, New Year Island"
  FT                   /collection_date="29-Mar-2023"
  ```

  `AC   ;` is a deliberate empty placeholder, not a missing accession. ENA requires exactly one `AC`
  block, and `AC   ;` and `AC * _<entry>` are not interchangeable.

- **Three facts are written differently in the flat file and the manifest.** Take all three from the
  `manifest` block, which uses Webin's vocabulary (Webin accepts only `PACBIO_SMRT` or `ILLUMINA` for
  `PLATFORM`):

  | Fact | Flat file | Manifest |
  |---|---|---|
  | assembler | `GetOrganelle v.1.7.7.0` (`CC Assembly Method`) | `GetOrganelle 1.7.7.0` |
  | platform | `Illumina` (`CC Sequencing Technology`) | `ILLUMINA` |
  | description | `… mitochondrion, complete genome` (`DE`) | `… mitochondrial genome` |

## `package_metadata.json`

| Field | Meaning |
|---|---|
| `schema_version` | `4` |
| `full_seqid` | the `ASSEMBLYNAME`; unique per assembly *and* annotation version |
| `og_id`, `assembly_prefix`, `annotation_version` | identity |
| `validation_study` | the `--ena_study` the run validated against. **Not a submission target**, see below |
| `mean_depth` | full-precision coverage; `manifest.COVERAGE` is its `%g` render |
| `mean_depth_source` | where `mean_depth` came from: `pipeline` (this run measured it), `database` (a stored value, e.g. a run with precomputed assemblies), `none`, or `unknown` (package built before this field existed) |
| `program`, `platform`, `scientific_name` | the values `manifest` renders, kept as this pipeline's record |
| `manifest` | the Webin genome-context keys this pipeline can fill, under Webin's names |
| `specimen` | the flat file `source` feature qualifiers |
| `sequence_length` | length in bases; matches the `ID` line |
| `sequence_sha256` | sha256 of the bases, uppercased, whitespace stripped |
| `normalised_circular_sha256` | rotation- and strand-invariant sequence identity |
| `flatfile_validation` | `status`, `reason`, `error_count`, `warning_count`, `webin_cli_version` |
| `package_digest` | changes when the submission content changes, and only then |

**Do not put `validation_study` in a manifest.** Sequence-context validation never resolves the
study, so nothing checks that it can receive data. The production launcher sets it to the
OceanOmics umbrella project `PRJEB110568`, which by definition cannot: a `STUDY` line copied from it
gets the submission rejected.

`package_digest` is a sha256 over the digests of the `.embl.gz`, `.chromosome_list.tsv.gz` and
`.tbl`, and nothing else. The gzipped files are written with `mtime=0` and an empty stored file name,
so their bytes, and the digest, are stable across reruns of the same assembly. `checksums.sha256` is
separate and covers every file in the directory except itself.

## More than one package per specimen

Every viable assembly and annotation version gets a package, and nothing deduplicates them, so
choosing between candidates belongs to the submission pipeline. Use `normalised_circular_sha256` for
equivalence: it is rotation- and strand-invariant, so the same circular genome assembled from a
different start base or on the other strand is recognised as the same sequence. Two candidates whose
digests match are interchangeable; two whose digests differ are a real disagreement, worth flagging
rather than resolving on recency.

## Locus tags

Nothing this pipeline writes carries a `/locus_tag`: not the `.tbl`, `.gbf`, `.embl.gz` or GFF. The
qualifier is **absent, not empty**, because an empty `/locus_tag=""` fails both table2asn and Webin.
Allocation and injection belong to the submission pipeline.

- Feature order in the flat file is stable and coordinate sorted, so numbering loci by walking the
  record gives the same answer on every rerun of the same assembly.
- table2asn reports `FATAL: NO_LOCUS_TAGS` for every record by design; this pipeline treats that code
  as advisory.

## BioSample

The package carries no BioSample, and the `SAMPLE` line belongs to the submission pipeline. One
thing to know when registering it: webin-cli resolves `SAMPLE` against ENA's own submission sample
service, which knows only samples registered through Webin, not the BioSamples mirror the ENA
browser serves. So `SAMEA` accessions resolve, and `SAMN` / `SAMD` (registered at NCBI or DDBJ) do
not. Verified against `ena-webin-cli 9.0.3`, `-context genome -validate -test`:

```text
SAMEA132129018 -> Submission(s) validated successfully.
SAMN40589646   -> Failed to initialise validator ... sample is null
```

The `specimen` block is the record to register the sample from: it is exactly what the flat file
carries, so the sample cannot contradict the assembly it belongs to.

## What `submission_ready` means

```sql
SELECT full_seqid, og_id, annotation, ena_study, webin_status, recorded_at
FROM ena_validation_attempts
WHERE submission_ready
ORDER BY og_id, full_seqid;
```

`submission_ready` means the flat file cleared every gate this pipeline runs, ending at the Webin
format check. It does **not** mean chosen, unembargoed, submittable or submitted: the study, sample
and runs are all still ahead of it. `ena_validation_latest` gives the most recent attempt per
`full_seqid`. Embargo comes from `sample.embargo_status` and is the submission pipeline's gate.

`ena_validation_attempts` rows are overwritten on rerun and carry no submission state, so they are
not a submission ledger. That is `ena_submissions`.

## Recording submissions

`ena_submissions` (`sql/015_ena_submissions.sql`) is the submission ledger. The submission pipeline
writes it and nothing else does. This pipeline never writes or reads it, deliberately: validation must
not depend on submission state, or a rerun of the validator could contradict the archive.

The key is `(full_seqid, webin_mode)`:

- **`full_seqid`**, not the assembly prefix, so a re-annotation is a new submittable record.
  `og_id`, `tech`, `seq_date`, `code` and `annotation` are `GENERATED ALWAYS` from `full_seqid`; do
  not write them.
- **`webin_mode`** is `production` or `test`, so a dry run against the test service cannot overwrite
  the production record of the same sequence.

Make the write idempotent with `ON CONFLICT ... DO UPDATE` and a predicate on `submission_status`,
so a rerun after a partial failure is safe:

```sql
INSERT INTO ena_submissions (
    full_seqid, webin_mode, submission_status, submitted_at, submitted_by,
    receipt_path, receipt_sha256, webin_cli_version,
    ena_study_accession, ena_analysis_accession,
    biosample_accession, biosample_source, locus_tag_prefix, run_accessions
) VALUES (...)
ON CONFLICT (full_seqid, webin_mode) DO UPDATE SET
    submission_status = EXCLUDED.submission_status,
    ena_analysis_accession = COALESCE(EXCLUDED.ena_analysis_accession,
                                      ena_submissions.ena_analysis_accession),
    updated_at = CURRENT_TIMESTAMP
WHERE ena_submissions.submission_status <> 'ACCESSION_ASSIGNED';
```

| Column | Holds |
|---|---|
| `submission_status` | `NOT_SUBMITTED` / `SUBMITTED` / `ACCESSION_ASSIGNED` / `FAILED`. Constraints enforce the evidence: `SUBMITTED` and `ACCESSION_ASSIGNED` need `submitted_at`, `ACCESSION_ASSIGNED` needs `ena_analysis_accession`, `FAILED` needs `error_message` |
| `submitted_at`, `submitted_by`, `webin_cli_version` | when, by whom, with which client |
| `receipt_path`, `receipt_sha256` | where the Webin receipt XML lives and its digest |
| `ena_study_accession` | `PRJEB…`, the BioProject submitted into |
| `ena_analysis_accession` | `ERZ…`, returned by webin-cli at submit time |
| `ena_assembly_accession`, `ena_sequence_accession` | `GCA_…` and the sequence accession, assigned later; filled on a second pass |
| `ena_sample_accession` | `ERS…`, only when a Webin sample was registered |
| `biosample_accession`, `biosample_source` | the BioSample the submission carried, and where it came from |
| `locus_tag_prefix` | the prefix registered with the study |
| `run_accessions` | the `ERR`/`SRR`/`DRR` runs the assembly came from |

The `ena_submission_status` view joins the latest validation attempt per `full_seqid` to the
production row of this table, so one query says what is validated and what has happened to it since.
Unsubmitted sequences read `NOT_SUBMITTED` rather than NULL.
