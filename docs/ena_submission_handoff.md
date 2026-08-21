# ENA submission handoff

For whoever maintains
[ENA-mito-genomes](https://github.com/Minderoo-OceanOmics-Centre-UWA/ENA-mito-genomes).

This describes the output the mitogenome pipeline publishes and how to submit from it.

## The seam

Specimens sit under embargo for an arbitrary period after their mitogenome is assembled, so
submission cannot be a stage of the assembly pipeline: the two events are months apart and the
second is not ours to trigger. The handover point is a directory on disk.

| | Mitogenomes (this repo) | ENA-mito-genomes |
|---|---|---|
| Builds the EMBL flatfile and chromosome list | yes | no |
| Format-validates the flatfile (`webin-cli -context sequence -validate`) | yes | no |
| Builds the submission manifest | no | yes |
| Allocates and injects locus tags | no | yes |
| Decides which assembly represents a specimen | no | yes |
| Applies the embargo gate | no | yes |
| Calls `webin-cli -context genome -submit` | no | yes |
| Records the resulting accessions | defines the table | writes it |

The mito pipeline only ever uses `-context sequence -validate` run
against the flatfile, and it validates the flatfile only, never submission metadata.
Genome-context validation cannot be run here at all: it resolves `SAMPLE` against ENA's
submission sample service before it checks anything else, and the sample is registered
downstream.

## Finding packages

```
s3://ocom-oceangenomes/analysed-data/mitogenomes/curated/<og_id>/<assembly_prefix>/
```

Each package describes itself; reading one needs no database join.

## What is in a package

Every file is named on `full_seqid`, so the whole package shares one stem.

| file | purpose |
|---|---|
| `<full_seqid>.embl.gz` | the flatfile. Webin `FLATFILE` input |
| `<full_seqid>.chromosome_list.tsv.gz` | Webin `CHROMOSOME_LIST` input |
| `<full_seqid>.package_metadata.json` | the handoff: everything this pipeline knows about the submission |
| `<full_seqid>.tbl` | NCBI feature table, kept for provenance |
| `<full_seqid>.fa` | sequence, collaborator handover |
| `<full_seqid>.gff` | annotation as the annotator produced it, collaborator handover |
| `<full_seqid>.genes.fa` | per-gene sequences (CDS, tRNA, rRNA), collaborator handover |
| `checksums.sha256` | sha256 of every other file in the directory |

The `.fa`, `.gff` and `.genes.fa` are **not** submission inputs. 

The chromosome list is a single row:

```
<full_seqid>	MT	Circular-Chromosome	Mitochondrion
```

## Building the manifest

The mito pipeline does not write a manifest. It cannot write a correct one: three of the
Webin genome-context keys are values only you hold.

- **`STUDY`** is the BioProject you register and submit into.
- **`SAMPLE`** is the BioSample you register for the specimen.
- **`RUN_REF`** names the raw-read submissions, which are yours.

So `package_metadata.json` is the whole handoff. Its `manifest` block holds the nine keys
this pipeline can fill, already under the names Webin uses, so building the manifest is
copying that block out as tab-separated lines and adding your own three. Injecting locus
tags into the flatfile or an accompanying list is the only edit needed to the package itself before `-submit`.

The `manifest` block from batch-01 OG51:

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

`FLATFILE` and `CHROMOSOME_LIST` are bare filenames, relative to the package directory.

A key whose value would not validate is left out rather than written wrong, so treat the
block as authoritative about what it contains and not as a fixed set of nine. In practice
all nine are present; a missing one means the pipeline had nothing trustworthy to put
there.

## Where each manifest value comes from

The manifest and the flatfile are independent webin-cli inputs. What the flatfile carries
does not license omitting anything from the manifest: every key above is required in the
manifest whatever the flatfile says, and the two are checked separately.

| manifest key | source |
|---|---|
| `STUDY` | yours, the registered BioProject |
| `SAMPLE` | yours, the registered BioSample |
| `RUN_REF` | yours, if ENA asks for run references |
| `ASSEMBLYNAME` | `manifest.ASSEMBLYNAME`, which is `full_seqid` |
| `ASSEMBLY_TYPE` | `manifest.ASSEMBLY_TYPE`, fixed at `clone or isolate` |
| `COVERAGE` | `manifest.COVERAGE`, the `%g` render of `mean_depth` |
| `PROGRAM` | `manifest.PROGRAM` |
| `PLATFORM` | `manifest.PLATFORM` |
| `MOLECULETYPE` | `manifest.MOLECULETYPE`, fixed at `genomic DNA` |
| `DESCRIPTION` | `manifest.DESCRIPTION`, built from `scientific_name` |
| `FLATFILE` | `manifest.FLATFILE` |
| `CHROMOSOME_LIST` | `manifest.CHROMOSOME_LIST` |

Notes that matter:

- **`full_seqid` is the `ASSEMBLYNAME`**, and it includes the annotation version, so a
  re-annotation is a new assembly rather than a collision with one you already sent. It is
  also the filename stem and the flatfile `ID` line identifier, so all three agree by
  construction.
- **You never need to re-derive a specimen fact.** The `specimen` block carries the
  qualifiers of the flatfile's `source` feature, read back out of the flatfile that was
  packaged, so the JSON and the flatfile cannot disagree. Neither the database nor the
  flatfile needs consulting for organism, isolate, tissue, place or date:

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

  A qualifier the specimen has no valid value for is **absent, not empty**: OG51 has no
  `/lat_lon`, so there is no `lat_lon` key. `organism` is normalised, so `Genus sp` reads `Genus sp.`.

- **The flatfile carries no study, sample, run or coverage.** The first three are not in
  the package either; they are yours. Coverage is `manifest.COVERAGE`. What the flatfile
  does carry:

  ```
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

  `AC   ;` is a deliberate empty placeholder, not a missing accession. ENA requires exactly
  one `AC` block, and `AC   ;` and `AC * _<entry>` are not interchangeable.

- **Three fields legitimately differ between the flatfile and the manifest**, so do not
  substitute one for the other. They are the same facts in two vocabularies:

  | fact | flatfile | manifest |
  |---|---|---|
  | assembler | `GetOrganelle v.1.7.7.0` (`CC Assembly Method`) | `GetOrganelle 1.7.7.0` |
  | platform | `Illumina` (`CC Sequencing Technology`) | `ILLUMINA` |
  | description | `… mitochondrion, complete genome` (`DE`) | `… mitochondrial genome` |

  Webin accepts only `PACBIO_SMRT` or `ILLUMINA` for `PLATFORM`. Take all three from the
  `manifest` block; they are already in Webin vocabulary there.

## `package_metadata.json`

| field | meaning |
|---|---|
| `schema_version` | `3` |
| `full_seqid` | the `ASSEMBLYNAME`; unique per assembly *and* annotation version |
| `og_id`, `assembly_prefix`, `annotation_version` | identity |
| `validation_study` | the study `--ena_study` named for the run, used only to make sequence-context validation execute. **Not a submission target**, see below |
| `mean_depth` | full-precision coverage; `manifest.COVERAGE` is its `%g` render |
| `program`, `platform`, `scientific_name` | the values `manifest` renders, kept as this pipeline's own record |
| `manifest` | the Webin genome-context keys this pipeline can fill, under Webin's names |
| `specimen` | the flatfile `source` feature qualifiers, so you need not parse the flatfile |
| `sequence_length` | length in bases, matches the `ID` line |
| `sequence_sha256` | sha256 of the bases, uppercased, whitespace stripped |
| `normalised_circular_sha256` | rotation- and strand-invariant sequence identity |
| `flatfile_validation` | nested: `status`, `reason`, `error_count`, `warning_count`, `webin_cli_version` |
| `package_digest` | changes when the ENA submission content changes, and only then |

**Do not put `validation_study` in a manifest.** It is whatever `--ena_study` was set to
when the run executed, and sequence-context validation never resolves it, so nothing
checks that it is submittable. It is `PRJEB110568`, the OceanOmics umbrella
project, which by definition cannot receive data: copying it into a `STUDY` line gets the
submission rejected. The field is recorded for provenance, not for reuse.

`package_digest` is a sha256 over the digest listing of the `.embl.gz`,
`.chromosome_list.tsv.gz` and `.tbl`. Those three files, and nothing else. The gzipped
files are written with `mtime=0` and an empty stored filename so their bytes are
reproducible and the digest is stable across reruns of the same assembly.
`checksums.sha256` is separate and covers every file in the directory except itself.

The webin-cli version that validated the flatfile is recorded in
`flatfile_validation.webin_cli_version`.

## BioSample

Registering the specimen is yours, end to end. This pipeline records no BioSample: the
package carries none, nothing here filters on one, and the `SAMPLE` line is yours to
write. What follows is what we learned about that step, because it is easy to get wrong.

webin-cli resolves `SAMPLE` against ENA's own submission sample service, which knows only
samples registered through Webin, not the EBI BioSamples mirror the ENA browser serves.
So:

- `SAMN` / `SAMD` means registered at another INSDC archive, and is **unusable at Webin**
- only `SAMEA` resolves

Verified against `ena-webin-cli 9.0.3`, `-context genome -validate -test`:

```
SAMEA132129018 -> Submission(s) validated successfully.
SAMN40589646   -> Failed to initialise validator ... sample is null
```

The `specimen` block in each package is the specimen record to register from: it is
exactly what was submitted in the flatfile, so a sample registered from it cannot
contradict the assembly it belongs to.

## More than one package per specimen

Every viable assembly and annotation version is published to the s3 bucket. Nothing deduplicates them, so choosing
between candidates is yours. Two things make it tractable:

- **Use `normalised_circular_sha256` for equivalence.** It is rotation- and strand-invariant, so a
  circular genome assembled starting at a different base, or on the other strand, is recognised as
  the same sequence rather than looking like a rival one. Two candidates whose digests match are
  interchangeable; two whose digests differ are a real disagreement and worth flagging rather than
  resolving silently on recency.

## No locus tags in anything we produce

Nothing this pipeline writes carries a `/locus_tag`: not the `.tbl`, not the `.gbf`, not the
`.embl.gz`, not the GFF. The qualifier is **absent, not empty**, because an empty `/locus_tag=""`
fails both `table2asn` and Webin. Allocation and injection are yours, end to end.

Two things that help:

- Feature order in the flatfile is stable and coordinate sorted. Numbering
  loci by walking the record gives the same answer on every rerun of the same assembly.
- `table2asn` reports `FATAL: NO_LOCUS_TAGS` for every record by design. This pipeline treats that
  code as advisory.

## What we vouch for

```sql
SELECT full_seqid, og_id, annotation, ena_study, webin_status, recorded_at
FROM ena_validation_attempts
WHERE submission_ready
ORDER BY og_id, full_seqid;
```

`submission_ready` means the flatfile cleared every gate this pipeline runs, ending at the Webin
format check. It does **not** mean chosen, unembargoed, submittable, or submitted: the study, the
sample and the runs are all still ahead of it. `ena_validation_latest` gives the most
recent attempt per `full_seqid` if you do not want to reason about attempt tokens.

Embargo comes from `sample.embargo_status` and is your gate; it should have exactly one owner.

## Writing back

There is now a table for you to close out: **`ena_submissions`** (`sql/015_ena_submissions.sql`).
It is written by you and nobody else. This pipeline never writes or reads it, deliberately:
validation must not depend on submission state, or a rerun of the validator could contradict the
archive.

Key is `(full_seqid, webin_mode)`:

- **`full_seqid`**, not the assembly prefix. It includes the annotation version, so a
  re-annotation is a new submittable thing rather than a collision with one you have already
  sent. `ena_validation_attempts` is keyed the same way, with the annotation version in its own
  `annotation` column. `og_id`, `tech`, `seq_date`, `code` and `annotation` on `ena_submissions`
  are `GENERATED ALWAYS` from `full_seqid`, so do not write them; they cannot drift.
- **`webin_mode`** is `production` or `test`. Keeping it in the key means a dry run against the
  test service cannot overwrite the production record of the same sequence.

**Make the write idempotent** with `ON CONFLICT (full_seqid, webin_mode) DO UPDATE` and a
predicate on `submission_status`. That is what makes a rerun after a partial failure safe:

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

| column | what it holds |
|---|---|
| `submission_status` | `NOT_SUBMITTED` / `SUBMITTED` / `ACCESSION_ASSIGNED` / `FAILED`. Constraints enforce the evidence: `SUBMITTED` and `ACCESSION_ASSIGNED` need `submitted_at`, `ACCESSION_ASSIGNED` needs `ena_analysis_accession`, `FAILED` needs `error_message` |
| `submitted_at`, `submitted_by`, `webin_cli_version` | when, by whom, with which client |
| `receipt_path`, `receipt_sha256` | where the Webin receipt XML lives and its digest, so the row can be traced back to the evidence on disk |
| `ena_study_accession` | `PRJEB…`, the BioProject you registered and submitted into |
| `ena_analysis_accession` | `ERZ…`, what webin-cli returns at submit time |
| `ena_assembly_accession`, `ena_sequence_accession` | `GCA_…` and `OU/LR…`, assigned later; fill them in on a second pass |
| `ena_sample_accession` | `ERS…`, only when a Webin sample was actually registered |
| `biosample_accession`, `biosample_source` | the BioSample **you registered or reused**, and where it came from. A `SAMN` here is a record that the specimen was NCBI-registered only, and so was not resolvable at Webin |
| `locus_tag_prefix` | the prefix registered at project creation |
| `run_accessions` | `ERR`/`SRR`/`DRR` runs the assembly came from |

`ena_submission_status` joins the latest validation attempt per `full_seqid` to the production row
of this table, so one query answers what is validated and what has happened to it since.
Un-submitted sequences read `NOT_SUBMITTED` rather than NULL.

Two things worth knowing about the current submitter before wiring it up:

- `07_push_submission_data_to_db.py` keys its write on `og_id/tech/seq_date/code`, with no
  annotation version, so two annotations of one assembly collide there. `ena_submissions` does
  not have that problem, but the value it writes has to be looked up by `full_seqid`.
- `04_create_biosample.py` skips Webin sample registration whenever `ncbi_biosample_id` is set,
  which is exactly the `SAMN` value webin-cli cannot resolve. Until that changes,
  `ena_sample_accession` stays NULL and `biosample_source` reads `sample.ncbi_biosample_id`,
  and the ledger will make that visible rather than leaving it implicit.

`ena_validation_attempts` rows are still overwritten freely on rerun and carry no submission
status, so do not treat them as a ledger; that is what `ena_submissions` is for.

`mitogenome_data.genbank_accession` is available if anything downstream reads it, but note the
column is named for GenBank and an ENA ERZ is not a GenBank accession.

## Questions for us

1. Does anything downstream read `mitogenome_data.genbank_accession`, or can the ERZ live only in
   your own tables? I think that `mitogenome_data.genbank_accession` can remain as is for now but
   all new ENA upload values live in our new table.
2. `package_metadata.json` is meant to be enough to build a manifest from, once you add
   `STUDY`, `SAMPLE` and `RUN_REF`. If a field is still missing, tell me: it is cheaper to
   add it to the metadata in the mito pipeline than for you to reconstruct it.
