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
| Builds the EMBL flatfile, allocates locus tags, builds the Webin package | yes | no |
| Validates the package locally and against Webin `-validate` | yes | no |
| Decides *which* assembly represents a specimen in a technology | yes | no |
| Applies the embargo gate | no | yes |
| Calls `webin-cli -submit` | no | yes |
| Records the resulting accessions | no | yes |

The seam between them is one database view.

## The contract

**We guarantee**, for every row in `ena_submission_queue`:

- an immutable package on disk at `package_path`, containing the gzipped EMBL flatfile, the
  gzipped chromosome list, a complete Webin `genome`-context manifest, and `checksums.sha256`
- locus tags already allocated, persisted in `ena_locus_registry`, and present in the flatfile
- the package has passed local validation and its manifest names a BioSample that Webin can resolve
- `full_seqid` is the `ASSEMBLYNAME` and is unique per assembly *and annotation version*, so a
  re-annotation is a new assembly rather than a collision

The package directory also holds `<full_seqid>.fa` and a locus-tagged `<full_seqid>.gff`. Those are
the collaborator handover pair, not submission inputs: the manifest never names them, and they are
excluded from `package_digest` so they cannot make an unchanged submission look changed.

**You guarantee**:

- you submit only what the queue lists, unmodified
- you close the row out afterwards (see [Writing back](#writing-back))

## Reading the queue

```sql
SELECT og_id, tech, full_seqid, ena_study_accession, biosample_accession,
       package_path, platform, mean_depth
FROM ena_submission_queue
WHERE embargo_status = 'Release'
ORDER BY og_id, tech;
```

The view already filters to `selection_status = 'SELECTED'`, `archive_status = 'NOT_SUBMITTED'`,
`package_status = 'READY'` and `local_validation_status = 'PASS'`. It exposes `embargo_status` but
deliberately does **not** filter on it, because embargo is your gate and should have exactly one
owner. Everything else is already decided.

Defined in `sql/008_ena_submission_queue.sql`.

### This replaces the samplesheet dedup

`workflow/00_create_mito_samplesheet.py:88-100` currently picks a candidate with:

```sql
ROW_NUMBER() OVER (PARTITION BY og_id, LOWER(tech)
                   ORDER BY seq_date_dt DESC NULLS LAST, seq_date DESC, code DESC)
```

That is the newest row in `mitogenome_data`. Our selection
(`bin/select_ena_submission.py`) instead picks among candidates that actually passed validation,
treats two assemblies as equivalent only when their `normalised_circular_sha256` matches (a
rotation- and strand-invariant digest, so a circular genome that starts at a different base is
recognised as the same sequence), and refuses to choose when they genuinely differ, marking the
specimen `MANUAL_REVIEW_REQUIRED`.

The two disagree whenever the newest assembly is not the one that passed. Today that means you can
submit a mitogenome nobody validated, and `MANUAL_REVIEW_REQUIRED` specimens are submitted anyway
because nothing tells you they are contested. Rows needing review never appear in the view, which is
the point of it.

Keep the `embargo_status = 'Release'` and `genbank_accession IS NULL` intent from the existing
query; drop the `ROW_NUMBER` dedup entirely.

## Writing back

Replace the update to `mitogenome_data.genbank_accession` in
`workflow/07_push_submission_data_to_db.py`:

```sql
UPDATE ena_submission_selections
SET archive_status        = 'SUBMITTED',
    ena_analysis_accession = %(erz)s,
    submitted_at           = CURRENT_TIMESTAMP,
    updated_at             = CURRENT_TIMESTAMP
WHERE ena_study_accession = %(study)s
  AND og_id               = %(og_id)s
  AND archive_status      = 'NOT_SUBMITTED';
```

and, once ENA assigns the assembly accession:

```sql
UPDATE ena_submission_selections
SET archive_status         = 'ACCESSION_ASSIGNED',
    ena_assembly_accession = %(gca)s,
    updated_at             = CURRENT_TIMESTAMP
WHERE ena_study_accession = %(study)s AND og_id = %(og_id)s;
```

The `archive_status = 'NOT_SUBMITTED'` predicate makes the write idempotent, so your existing
`UPDATED` / `SKIPPED` / `CONFLICT_SKIPPED` / `ROW_NOT_FOUND` outcomes carry over unchanged.

This write is load-bearing beyond bookkeeping. `select_ena_submission.py` refuses to move a
selection whose `archive_status` is anything other than `NOT_SUBMITTED`, so until you set it, a
re-run of selection can silently repoint a specimen you have already submitted.

Keep writing `mitogenome_data.genbank_accession` too if anything downstream reads it, but note that
the column is named for GenBank and an ENA ERZ is not a GenBank accession.

## Steps that become lookups

With the per-technology studies live (below), three of the eight steps are no longer needed:

- **03, project registration.** The studies exist. Read them from `assets/ena/accessions.tsv`.
- **04, BioSample registration.** The manifest BioSample is already resolved and recorded in the
  queue. See [BioSample](#biosample) for the part that is genuinely unsolved.
- **05, manifest building and locus-tag injection.** The package is built. Read it, do not rebuild
  it.

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

|  | this repo | ENA-mito-genomes |
|---|---|---|
| prefix | per technology, `OGMTHIFI` | the OG ID, `OG910` |
| rendered tag | `OGMTHIFI_000910001` | `OG910_H00001` |
| allocated by | `bin/allocate_ena_locus_tags.py`, DB registry under `pg_advisory_xact_lock` | in-memory, order of first `/gene=` appearance (`05:162-198`) |
| injected into | the NCBI `.tbl`, before `table2asn` | the EMBL flatfile, after annotation |
| stable across reruns | yes, persisted | no, positional |

`sql/006_ena_tech_aware_locus_tags.sql` constrains `ena_candidate_loci.locus_tag` to
`^[A-Z][A-Z0-9]{2,11}_[0-9]{9}$`, which `OG910_H00001` does not satisfy. Injecting into the EMBL
also means `table2asn` never sees the tags and cannot validate them.

Since tags are allocated long before embargo lifts, read them rather than deriving them.

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
`package_status = 'BLOCKED_NCBI_ONLY_BIOSAMPLE'`, which keeps those specimens out of the queue
rather than letting them fail at submission time. `bin/audit_ena_readiness.py --check-ena-mirror`
queries the EBI browser API per specimen if you want the current picture.

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
   `ena_submission_selections`?
2. Do you want the queue to expose anything else, so you never have to join back to
   `ena_candidate_packages` yourself?
