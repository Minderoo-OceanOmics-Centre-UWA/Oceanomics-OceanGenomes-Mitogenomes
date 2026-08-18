# ENA submission handoff

For whoever maintains
[ENA-mito-genomes](https://github.com/Minderoo-OceanOmics-Centre-UWA/ENA-mito-genomes).

## Why there are two pipelines

Specimens sit under embargo for an arbitrary period after their mitogenome is assembled, so
submission cannot be a stage of the assembly pipeline: the two events are months apart and the
second is not ours to trigger. The split is deliberate and stays.

What follows from that is a division of ownership:

| | Mitogenomes (this repo) | ENA-mito-genomes |
|---|---|---|
| Builds the EMBL flatfile and the Webin package | yes | no |
| Allocates locus tags and injects them into the flatfile | no | yes |
| Format-validates the flatfile with Webin `-validate -context sequence` | yes | no |
| Validates the package with Webin `-validate -context genome` | no | yes |
| Decides *which* assembly represents a specimen in a technology | no | yes |
| Applies the embargo gate | no | yes |
| Calls `webin-cli -submit` | no | yes |
| Records the resulting accessions | no | yes |

The seam between them is the published package directory on disk, plus one
boolean in the database.

> **Changed 2026-08.** Selection used to live here, behind an
> `ena_submission_queue` view backed by `ena_candidate_packages` and
> `ena_submission_selections`. `sql/012_drop_ena_selection_layer.sql` drops all
> three: choosing among candidates is a submission-time decision and belongs
> with the embargo gate, in one place. What this repo now asserts is narrower and
> per-candidate: this flatfile is well formed. The sections below describe the
> replacement; the git history has the old contract if you need it.

## The contract

**We guarantee**, for every published package directory whose validation row is
`submission_ready`:

- an immutable package on disk, containing the gzipped EMBL flatfile, the
  gzipped chromosome list, a Webin `genome`-context manifest, `checksums.sha256`, and
  `<full_seqid>.package_metadata.json`
- a flatfile carrying **no** `/locus_tag` qualifier on any feature: tag allocation and injection
  are yours, and the qualifier is absent rather than empty because an empty one fails validation
- a flatfile that passed `ena-webin-cli -context sequence -validate`, which is what
  `submission_ready` records
- `full_seqid` is the `ASSEMBLYNAME` and is unique per assembly *and annotation version*, so a
  re-annotation is a new assembly rather than a collision

**We no longer guarantee**, and this is the part to read twice:

- **a resolvable BioSample.** The old queue filtered to `package_status = 'READY'`, so a package
  whose accession was missing or NCBI-only never reached you. Nothing filters now. Read
  `biosample_accession` from `package_metadata.json` and check it yourself: absent means
  unregistered, `SAMN`/`SAMD` means registered at another INSDC archive and unusable at Webin, and
  only `SAMEA` resolves. See [BioSample](#biosample) for why. A blocked package's manifest starts
  with `# BLOCKED:` and carries no `SAMPLE` line, so it is also self-describing on disk.
- **one package per specimen per technology.** Every viable assembly and annotation version is
  published. Deduplicating them is now yours; see [Choosing among candidates](#choosing-among-candidates).

The package directory also holds `<full_seqid>.fa`, `<full_seqid>.gff` and `<full_seqid>.tbl`. The
FASTA and GFF are the collaborator handover pair, not submission inputs: the manifest never names
them, and they are excluded from `package_digest` so they cannot make an unchanged submission look
changed. The GFF is the annotator's own output, untouched; it carries no `locus_tag=` attribute for
the same reason the flatfile does not.

**You guarantee**:

- you submit packages unmodified, apart from injecting locus tags
- you close the row out afterwards (see [Writing back](#writing-back))

## Finding packages

Packages are published to:

```
<outdir>/mitogenomes/<og_id>/<assembly_prefix>/ena/package/
```

so a glob over `*/*/ena/package/*.package_metadata.json` enumerates every candidate. Each metadata
file is the package's own description and needs no join:

| field | use |
|---|---|
| `full_seqid` | the Webin `ASSEMBLYNAME`; unique per assembly *and* annotation version |
| `og_id`, `assembly_prefix`, `annotation_version` | identity |
| `study` | the ENA child study for this technology |
| `biosample_accession` | check it yourself; see the contract above |
| `platform`, `program`, `mean_depth`, `scientific_name`, `run_accessions` | manifest inputs, already resolved |
| `sequence_sha256`, `normalised_circular_sha256` | equivalence; see below |
| `package_digest` | changes when the ENA submission content changes, and only then |
| `flatfile_validation` | status, reason, error/warning counts, webin-cli version |

To restrict to candidates this pipeline vouches for, join on the validation table:

```sql
SELECT assembly_prefix, og_id, ena_study, webin_status, recorded_at
FROM ena_validation_attempts
WHERE submission_ready
ORDER BY og_id, assembly_prefix;
```

`submission_ready` means the flatfile cleared every gate this pipeline runs, ending at the Webin
format check. It does not mean chosen, unembargoed, or submitted. `ena_validation_latest` gives the
most recent attempt per `assembly_prefix` if you do not want to reason about attempt tokens.
Embargo comes from `sample.embargo_status`, and is your gate: it should have exactly one owner.

### Choosing among candidates

`workflow/00_create_mito_samplesheet.py:88-100` currently picks a candidate with:

```sql
ROW_NUMBER() OVER (PARTITION BY og_id, LOWER(tech)
                   ORDER BY seq_date_dt DESC NULLS LAST, seq_date DESC, code DESC)
```

That is the newest row in `mitogenome_data`, which disagrees with validation whenever the newest
assembly is not the one that passed. Selecting on recency alone means you can submit a mitogenome
nobody validated.

The approach this repo used before selection moved to you, worth reusing rather than rederiving:

- consider only candidates whose validation row is `submission_ready`
- select within a technology, never across: a specimen can legitimately be submitted once as HiFi,
  once as Hi-C and once as Illumina, because those are separate assemblies from separate data
- treat two candidates as the same sequence only when their `normalised_circular_sha256` matches.
  That digest is rotation- and strand-invariant, so a circular genome assembled starting at a
  different base is recognised as the same sequence rather than looking like a rival one
- among equivalent candidates, take the newest sequencing date
- when candidates genuinely differ, do not pick: flag the specimen for review. Silently taking the
  newest is how a contested specimen gets submitted without anyone noticing

Keep the `embargo_status = 'Release'` and `genbank_accession IS NULL` intent from the existing
query; drop the `ROW_NUMBER` dedup entirely.

## Writing back

There is no longer a table here for you to close out: `ena_submission_selections` went with the
selection layer, so submission state is yours to keep, in whatever schema suits you. Two things
worth carrying over from the design that was there:

- key your own record on `full_seqid`, not `assembly_prefix`. `full_seqid` includes the annotation
  version, so a re-annotation is a new submittable thing rather than a collision with one you have
  already sent.
- make the write idempotent with a predicate on your own status column, the way the
  `archive_status = 'NOT_SUBMITTED'` predicate did. That is what makes a rerun after a partial
  failure safe, and it preserves your existing `UPDATED` / `SKIPPED` / `CONFLICT_SKIPPED` /
  `ROW_NOT_FOUND` outcomes unchanged.

Nothing in this repo reads that state back. `ena_validation_attempts` rows are overwritten freely on
rerun and carry no submission status, so do not treat them as a submission ledger.

Keep writing `mitogenome_data.genbank_accession` too if anything downstream reads it, but note that
the column is named for GenBank and an ENA ERZ is not a GenBank accession.

## Steps that become lookups

With the per-technology studies live (below), three of the eight steps are no longer needed:

- **03, project registration.** The studies exist. Read them from `assets/ena/accessions.tsv`.
- **04, BioSample registration.** The manifest BioSample is already resolved and recorded in
  `package_metadata.json`, though verifying it resolves at Webin is now yours. See
  [BioSample](#biosample) for the part that is genuinely unsolved.
- **05, manifest building.** The package and its manifest are built. Read them, do not rebuild
  them. Locus-tag injection is the one part of step 05 that stays yours.

Steps 01, 02 and 06 survive, as does the duplicate-alias reuse logic (see
[Worth keeping](#worth-keeping)).

## The project structure, which needs settling

`assets/ena/prod_add_20260806T133137.xml` in this repo is a real ENA receipt, `success="true"`,
dated 2026-08-06:

- `PRJEB110568` is an **umbrella** project (parent `PRJNA1046164`) and cannot receive data
- three child studies were created, one per sequencing technology, held to 2028-08-06:

  | tech | study | locus tag prefix |
  |---|---|---|
  | hifi | `PRJEB123419` | `OGMTHIFI` |
  | hic  | `PRJEB123420` | `OGMTHIC` |
  | ilmn | `PRJEB123421` | `OGMTILMN` |

`workflow/03_create_og_bioproject.py:204-229` registers **one PROJECT per OG** under that same
umbrella, with a `LOCUS_TAG_PREFIX` equal to the OG ID (`03:95-104`). At catalogue scale that is
thousands of child studies and thousands of globally-unique ENA locus tag prefixes claimed, against
an umbrella that already has its three children, from a second codebase.

Your `receipts/` is empty and `nextflow-run.sh:17` points at `wwwdev`, so this looks un-run in
production. Please keep it that way until we have agreed the structure. This is the one item here
that cannot be resolved by changing code on either side alone.

## Locus tags

**Locus tags are yours, end to end.** This repo used to allocate them and inject them into the NCBI
`.tbl` before `table2asn`; it no longer does. Nothing it produces carries a `/locus_tag`: not the
`.tbl`, not the `.gbf`, not the `.embl.gz`, not the collaborator GFF. The qualifier is absent, not
empty, because an empty `/locus_tag=""` fails both `table2asn` and Webin.

What that means in practice:

- Inject tags into the flatfile after you read the package, before `webin-cli -submit`.
- The prefixes registered against the three child studies are in the table above and in
  `assets/ena/accessions.tsv`. They stay registered and they are the ones to use; this repo simply
  no longer resolves or applies them.
- Feature order in the flatfile is stable and coordinate sorted. `bin/process_files.py`
  (`sort_tbl_features`) re-sorts Emma's string-ordered output numerically, so numbering loci by
  walking the record gives the same answer on every rerun of the same assembly.
- `table2asn` reports `FATAL: NO_LOCUS_TAGS` for every record by design. This repo treats that code
  as advisory (`bin/parse_table2asn_validation.py`); it is not a defect in the package.
- The old registry tables (`ena_locus_registry`, `ena_candidate_loci`) are dropped by
  `sql/010_drop_ena_locus_tables.sql`, with their contents preserved in `*_archive` tables. Those
  archives are frozen history for records already submitted. Do not read them as a live source of
  tags.

## BioSample

`04_create_biosample.py:263-269` skips SAMPLE creation whenever `ncbi_biosample_id` is set, which is
the common case, and assumes the resulting `SAMN` works at ENA. It does not, and this is the one
genuinely open problem on both sides.

webin-cli resolves the manifest `SAMPLE` against ENA's own submission sample service, which knows
only samples registered through Webin, not the EBI BioSamples mirror the ENA browser serves.
Verified against `ena-webin-cli 9.0.3`, `-context genome -validate -test`:

```
SAMEA132129018 -> Submission(s) validated successfully.
SAMN40589646   -> Failed to initialise validator ... sample is null
```

See `sql/007_insdc_biosample_accessions.sql` for the full note. We widened the accession CHECK to
`SAM(EA|N|D)` so the catalogue can be recorded, and express the unusability as
`SAM(EA|N|D)` so the catalogue can be recorded. This repo used to express the unusability as
`package_status = 'BLOCKED_NCBI_ONLY_BIOSAMPLE'` and filter those specimens out of the queue; with
the queue gone, that check is yours, and the `# BLOCKED:` manifest is the on-disk signal.
`bin/audit_ena_readiness.py --check-ena-mirror` queries the EBI browser API per specimen if you want
the current picture.

## Manifest details worth not re-deriving

Recorded so the differences read as reasons rather than preferences. Each is a Webin failure that
step 05 would reintroduce if the manifest were rebuilt:

| | `bin/ena_package.py` | `05_make_assembly_manifests.py` |
|---|---|---|
| `PLATFORM` | `PACBIO_SMRT` / `ILLUMINA` | `PacBio` / `Illumina` (`:220-226`), not Webin's vocabulary |
| organism | `normalise_open_nomenclature` rewrites `Genus sp` to `Genus sp.` | verbatim (`03:171-198`); webin rejects with "Organism is not Submittable" |
| chromosome list | `MT / Circular-Chromosome / Mitochondrion` (`:148`) | `MIT / chromosome / Mitochondrion` (`:257-261`); never declares circularity |
| manifest keys | includes `MOLECULETYPE`, `DESCRIPTION`, `RUN_REF` (`:153`) | omits all three |
| `ASSEMBLYNAME` | 5-token `full_seqid` | 4-token prefix (`:431`); re-annotation collides |
| lat/lon | from `bin/build_source_modifiers.py` | declares `UNITS="DD"` (`04:118-119`) over DMS text such as `114o 15.765`; ERC000011 rejects it |

Two more, unrelated to the manifest:

- `06:152-164` calls `-submit` directly, with no `-validate` pass first. Worth adding regardless of
  anything here.
- `02:178-184` exits 0 on validation failure and drops the failing pairs, so a run where nothing
  validated looks like a successful run. The last recorded run was 0/1 passing and steps 03-07 never
  executed.

Unrelated to any of this: `README.md:16-19` has an unresolved merge conflict with no closing marker.

## Worth keeping

Things that repo does that this one does not, and that should survive the refactor:

- the webin-cli container pinned by sha256 digest (`nextflow.config:19`)
- reading the receipt from `genome/*/submit/receipt.xml` rather than `submission.xml`
  (`06:70-87`) and preserving the webin output tree *before* checking the return code
  (`06:185-188`), so failures stay diagnosable
- parsing ENA's duplicate-alias error to reuse an existing accession, which is what makes re-runs
  idempotent. Currently implemented three times (`03:106-147`, `04:204-234`, `05:288-332`); one
  copy would do
- the conflict-safe writeback semantics in step 07
- the embargo gate itself

## Questions for us

1. Does anything downstream read `mitogenome_data.genbank_accession`, or can the ERZ live only in
   your own tables?
2. Is `package_metadata.json` plus `submission_ready` enough to build a manifest from, or is there
   a field you would otherwise have to rederive? It is cheaper to add it to the metadata here than
   for you to reconstruct it.
