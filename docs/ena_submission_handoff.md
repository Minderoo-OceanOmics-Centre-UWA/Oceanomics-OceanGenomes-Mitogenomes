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
| Builds the EMBL flatfile, chromosome list and Webin manifest | yes | no |
| Format-validates the flatfile (`webin-cli -context sequence -validate`) | yes | no |
| Allocates and injects locus tags | no | yes |
| Decides which assembly represents a specimen | no | yes |
| Applies the embargo gate | no | yes |
| Calls `webin-cli -context genome -submit` | no | yes |
| Records the resulting accessions | defines the table | writes it |

## Finding packages

```
<outdir>/mitogenomes/<og_id>/<assembly_prefix>/ena/package/
```

A glob over `*/*/ena/package/*.package_metadata.json` enumerates every candidate. Each package
describes itself; no database join is needed to build a manifest from it.

## What is in a package

Every file is named on `full_seqid`, so the whole package shares one stem.

| file | purpose |
|---|---|
| `<full_seqid>.embl.gz` | the flatfile. Webin `FLATFILE` input |
| `<full_seqid>.chromosome_list.tsv.gz` | Webin `CHROMOSOME_LIST` input |
| `<full_seqid>.manifest.txt` | the Webin `genome`-context manifest, already built |
| `<full_seqid>.package_metadata.json` | every value the manifest was built from |
| `<full_seqid>.tbl` | NCBI feature table, kept for provenance |
| `<full_seqid>.fa` | sequence, collaborator handover |
| `<full_seqid>.gff` | annotation as the annotator produced it, collaborator handover |
| `<full_seqid>.genes.fa` | per-gene sequences (CDS, tRNA, rRNA), collaborator handover |
| `checksums.sha256` | sha256 of every other file in the directory |

The `.fa`, `.gff` and `.genes.fa` are **not** submission inputs. The manifest never names them and
they are excluded from `package_digest`, so a change to them cannot make an unchanged submission
look changed.

The chromosome list is a single row:

```
<full_seqid>	MT	Circular-Chromosome	Mitochondrion
```

## The manifest is already built

Read it, do not rebuild it. Injecting locus tags into the flatfile is the only edit needed before
`-submit`.

```
STUDY	PRJEB123419
SAMPLE	SAMEA132129018
ASSEMBLYNAME	OG2124.hifi.260421.v323mitohifi.emma102
ASSEMBLY_TYPE	clone or isolate
COVERAGE	192.749
PROGRAM	MitoHiFi 3.2.3
PLATFORM	PACBIO_SMRT
MOLECULETYPE	genomic DNA
DESCRIPTION	Ateleopus japonicus mitochondrial genome
FLATFILE	OG2124.hifi.260421.v323mitohifi.emma102.embl.gz
CHROMOSOME_LIST	OG2124.hifi.260421.v323mitohifi.emma102.chromosome_list.tsv.gz
```

`RUN_REF` is emitted when the package carries run accessions. `run_accessions` is empty in the
packages we publish, so if ENA wants run references, supplying them is yours.

## Where each manifest value comes from

If you build your own manifest rather than reading ours, this is the mapping. The short version:
**take submission identifiers and the Webin-vocabulary strings from the JSON; the flatfile already
carries everything about the specimen itself.**

| manifest key | `package_metadata.json` field | in the flatfile? |
|---|---|---|
| `STUDY` | `study` | no |
| `SAMPLE` | `biosample_accession` (check it, see below) | no |
| `ASSEMBLYNAME` | `full_seqid` | yes, `ID` line and `AC * _` entry name |
| `ASSEMBLY_TYPE` | fixed, `clone or isolate` | no |
| `COVERAGE` | `mean_depth` | no |
| `PROGRAM` | `program` | only as free text, see the warning below |
| `PLATFORM` | `platform` | only as free text, see the warning below |
| `MOLECULETYPE` | fixed, `genomic DNA` | yes, `ID` line and `/mol_type` |
| `DESCRIPTION` | built from `scientific_name` | not verbatim, see the warning below |
| `FLATFILE` | `<full_seqid>.embl.gz` | n/a |
| `CHROMOSOME_LIST` | `<full_seqid>.chromosome_list.tsv.gz` | n/a |

Notes that matter:

- **`full_seqid` is the `ASSEMBLYNAME`**, and it includes the annotation version, so a
  re-annotation is a new assembly rather than a collision with one you already sent. It is also
  the filename stem and the flatfile `ID` line identifier, so all three agree by construction.
- **The flatfile already carries every specimen fact**, so none of it needs re-deriving: organism
  (`OS` and `/organism`, normalised so `Genus sp` reads `Genus sp.`), topology (`circular` on the
  `ID` line), molecule type, sequence length, and a full `source` feature with `/organelle`,
  `/isolate`, `/tissue_type`, `/geo_loc_name`, `/lat_lon` and `/collection_date`. `/lat_lon` and
  `/collection_date` are omitted rather than filled with a placeholder when the specimen has no
  valid value.

  ```
  ID   OG2124.hifi.260421.v323mitohifi.emma102; SV 1; circular; genomic DNA; STD; UNC; 16650 BP.
  XX
  AC   ;
  XX
  AC * _OG2124.hifi.260421.v323mitohifi.emma102
  XX
  DE   Ateleopus japonicus mitochondrion, complete genome
  XX
  OS   Ateleopus japonicus
  XX
  CC   ##Assembly-Data-START##
  CC   Assembly Method       :: MitoHifi v.3.2.3
  CC   Sequencing Technology :: PacBio HiFi
  CC   ##Assembly-Data-END##
  ...
  FT   source          1..16650
  FT                   /organism="Ateleopus japonicus"
  FT                   /organelle="mitochondrion"
  FT                   /mol_type="genomic DNA"
  FT                   /isolate="OG2124"
  FT                   /tissue_type="Gills"
  FT                   /geo_loc_name="Australia: EEZ"
  FT                   /lat_lon="21.57628 S 156.49365 E"
  FT                   /collection_date="25-Oct-2025"
  ```

  `AC   ;` is a deliberate empty placeholder, not a missing accession. ENA requires exactly one
  `AC` block, and `AC   ;` and `AC * _<entry>` are not interchangeable.

- **The flatfile carries no study, sample, run or coverage.** Those exist only in the JSON.
- **Do not substitute the flatfile's `CC ##Assembly-Data##` strings for `PROGRAM` and
  `PLATFORM`.** They describe the same facts in a different vocabulary: the flatfile says
  `MitoHifi v.3.2.3` and `PacBio HiFi`, the manifest needs `MitoHiFi 3.2.3` and `PACBIO_SMRT`.
  Webin only accepts `PACBIO_SMRT` or `ILLUMINA`. Likewise the flatfile `DE` line
  (`… mitochondrion, complete genome`) is not the manifest `DESCRIPTION`
  (`… mitochondrial genome`). Take all three from the JSON.

## `package_metadata.json`

| field | meaning |
|---|---|
| `schema_version` | `2` |
| `full_seqid` | the `ASSEMBLYNAME`; unique per assembly *and* annotation version |
| `og_id`, `assembly_prefix`, `annotation_version` | identity |
| `study` | the ENA child study recorded for this assembly's technology |
| `biosample_accession` | may be `null`; check it yourself, see below |
| `mean_depth` | manifest `COVERAGE` |
| `program`, `platform`, `scientific_name` | manifest inputs in Webin vocabulary |
| `run_accessions` | empty; run references are yours to add |
| `sequence_length` | length in bases, matches the `ID` line |
| `sequence_sha256` | sha256 of the bases, uppercased, whitespace stripped |
| `normalised_circular_sha256` | rotation- and strand-invariant sequence identity |
| `flatfile_validation` | nested: `status`, `reason`, `error_count`, `warning_count`, `webin_cli_version` |
| `package_digest` | changes when the ENA submission content changes, and only then |

`package_digest` is a sha256 over the digest listing of the `.embl.gz`, `.chromosome_list.tsv.gz`,
`.manifest.txt` and `.tbl`. The gzipped files are written with `mtime=0` and an empty stored
filename so their bytes are reproducible and the digest is stable across reruns of the same
assembly. `checksums.sha256` is separate and covers every file in the directory except itself.

The webin-cli version that validated the flatfile is recorded in
`flatfile_validation.webin_cli_version`.

## BioSample

**Read `biosample_accession` and check it before you submit.** Nothing upstream filters on it.

webin-cli resolves the manifest `SAMPLE` against ENA's own submission sample service, which knows
only samples registered through Webin, not the EBI BioSamples mirror the ENA browser serves. So:

- absent or `null` means unregistered
- `SAMN` / `SAMD` means registered at another INSDC archive and **unusable at Webin**
- only `SAMEA` resolves

Verified against `ena-webin-cli 9.0.3`, `-context genome -validate -test`:

```
SAMEA132129018 -> Submission(s) validated successfully.
SAMN40589646   -> Failed to initialise validator ... sample is null
```

A package that fails this check is self-describing on disk. Its manifest starts with a
`# BLOCKED:` line carrying the reason, and keeps only four keys, with no `SAMPLE`:

```
# BLOCKED: Invalid or missing ENA BioSample accession: ''
STUDY	PRJEB123419
ASSEMBLYNAME	OG2124.hifi.260421.v323mitohifi.emma102
FLATFILE	OG2124.hifi.260421.v323mitohifi.emma102.embl.gz
CHROMOSOME_LIST	OG2124.hifi.260421.v323mitohifi.emma102.chromosome_list.tsv.gz
```

**Expect this to be the common case.** As of 2026-08-20, all 219 published packages under
`/scratch/pawsey1348/tpeirce` are blocked: `biosample_accession` is `null` in 202 of them and
`SAMN…` in the rest, and none is `SAMEA`. Until specimens are registered through Webin, no
package is submittable. `bin/audit_ena_readiness.py --check-ena-mirror` queries the EBI browser
API per specimen if you want the current picture.

## More than one package per specimen

Every viable assembly and annotation version is published. Nothing deduplicates them, so choosing
between candidates is yours. Two things make it tractable:

- **Select within a technology, never across.** A specimen can legitimately be submitted once as
  HiFi, once as Hi-C and once as Illumina; those are separate assemblies from separate data.
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

- Feature order in the flatfile is stable and coordinate sorted. `bin/process_files.py`
  (`sort_tbl_features`) re-sorts the annotator's string-ordered output numerically, so numbering
  loci by walking the record gives the same answer on every rerun of the same assembly.
- `table2asn` reports `FATAL: NO_LOCUS_TAGS` for every record by design. This pipeline treats that
  code as advisory (`bin/parse_table2asn_validation.py`); it is not a defect in the package.

## What we vouch for

```sql
SELECT full_seqid, og_id, annotation, ena_study, webin_status, recorded_at
FROM ena_validation_attempts
WHERE submission_ready
ORDER BY og_id, full_seqid;
```

`submission_ready` means the flatfile cleared every gate this pipeline runs, ending at the Webin
format check. It does **not** mean chosen, unembargoed, submittable, or submitted: a
`submission_ready` package can still be BioSample-blocked. `ena_validation_latest` gives the most
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
| `ena_study_accession` | `PRJEB…`, the per-technology child study it actually went to |
| `ena_analysis_accession` | `ERZ…`, what webin-cli returns at submit time |
| `ena_assembly_accession`, `ena_sequence_accession` | `GCA_…` and `OU/LR…`, assigned later; fill them in on a second pass |
| `ena_sample_accession` | `ERS…`, only when a Webin sample was actually registered |
| `biosample_accession`, `biosample_source` | the BioSample the **manifest carried**, and where it came from. A `SAMN` here is a record that the specimen was NCBI-registered only |
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
  `ena_sample_accession` stays NULL and `biosample_source` reads `sample.ncbi_biosample_id` —
  and the ledger will make that visible rather than leaving it implicit.

`ena_validation_attempts` rows are still overwritten freely on rerun and carry no submission
status, so do not treat them as a ledger; that is what `ena_submissions` is for.

`mitogenome_data.genbank_accession` is available if anything downstream reads it, but note the
column is named for GenBank and an ENA ERZ is not a GenBank accession.

## Questions for us

1. Does anything downstream read `mitogenome_data.genbank_accession`, or can the ERZ live only in
   your own tables?
2. Is `package_metadata.json` enough to build a manifest from, or is there a field you would
   otherwise have to rederive? It is cheaper to add it to the metadata here than for you to
   reconstruct it.
