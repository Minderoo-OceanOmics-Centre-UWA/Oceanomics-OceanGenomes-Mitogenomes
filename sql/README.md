# Database migrations

Numbered, forward-only migrations for the OceanOmics PostgreSQL schema this pipeline reads and
writes. Apply them with `bin/apply_ena_migrations.py`, never by hand; see
[Database migrations](../docs/usage.md#database-migrations) for how the runner works. Each file's
header comment has the full background.

| Migration | What it does |
|---|---|
| `001_create_ena_validation_attempts` | Creates `ena_validation_attempts`, one record per ENA validation |
| `002_ena_validation_attempts_single_row_per_attempt` | One row per attempt key, so reruns overwrite instead of appending |
| `003_mitogenome_data_uniform_depth` | Adds the uniform remap depth (`mean_depth`, `depth_method`) to `mitogenome_data` and relabels legacy rows |
| `004_ena_candidate_packages` | Package registry and locus-tag tables (later retired by 010 to 013) |
| `005_ena_validation_attempts_genome_context` | Splits test/production Webin results (columns later retired by 011 and 012) |
| `006_ena_tech_aware_locus_tags` | Per-technology locus-tag prefixes (later retired by 010) |
| `007_insdc_biosample_accessions` | Records any INSDC BioSample accession on `ena_specimen_accessions` and `ena_candidate_packages` |
| `008_ena_submission_queue` | Read-only queue view for the submitter (later retired by 012) |
| `009_ena_locus_registry_canonical_order` | Renumbers locus serials into coordinate order (later retired by 010) |
| `010_drop_ena_locus_tables` | Retires the locus-tag registry; archives it first |
| `011_drop_local_package_validation` | Retires the local package validation gate |
| `012_drop_ena_selection_layer` | Retires the in-repo selection layer; **not reversible** |
| `013_drop_ena_candidate_runs` | Retires `ena_candidate_runs` |
| `014_ena_validation_attempts_full_seqid` | Re-keys `ena_validation_attempts` on `full_seqid`; table rebuild |
| `015_ena_submissions` | Adds the `ena_submissions` ledger and `ena_submission_status` view |
| `016_mitogenome_data_og_num_generated` | Makes `mitogenome_data.og_num` a stored generated column |
| `017_lca_content_addressed_rows` | Records the content-addressed keys on `lca` / `lca_raw_results` |
| `018_mitogenome_data_og_num_first` | Moves `og_num` back to column 1 of `mitogenome_data`; table rebuild |
| `019_ena_validation_attempts_recompute_submission_ready` | Restores `submission_ready` on rows a replay of 005 wrongly cleared |
| `020_ena_validation_attempts_og_num` | Adds generated `og_num` to `ena_validation_attempts`; table rebuild |
| `021_mitogenome_data_trna_advisory` | Records the tRNAs the annotation QC gate waived (`trna_advisory`) |
| `022_lca_widen_taxon_id_and_confidence` | Widens LCA columns that were too narrow for real BLAST output |
| `023_mitogenome_data_order_variant` | Records why a non-canonical gene order was accepted, or how far off a rejected one was |
| `024_lca_validation_validated_rank` | Records the rank at which a species ID was validated |
| `025_mitogenome_data_order_variant_taxon_check` | Records whether a second taxonomy source agreed with a gene-order rule |
| `026_mitogenome_data_annotation_integrity` | Additive columns for the final annotation-integrity verdict and conversion audit |

The runner replays the whole list against a fresh database, so the early ENA migrations still run
even though 010 to 013 later remove much of what they create.

## Migrations that need care

**010** retires `ena_locus_registry` and `ena_candidate_loci` now that locus tags are assigned by
the downstream submission pipeline. It copies both into `*_archive` tables first: the tags behind
already-submitted records cannot be rebuilt from the flat files.

**012** retires the selection layer for the same reason: choosing and submitting a package belongs
to the downstream pipeline, so `ena_candidate_packages`, `ena_submission_selections` and the
`ena_submission_queue` view are dropped, along with the package, production-Webin, checksum and
run-provenance columns on `ena_validation_attempts`. That leaves `submission_ready` meaning what
this pipeline can actually attest: the flat file passed every gate. The drop is not reversible from
the repository, so dump those three relations first if you want them.

**014** re-keys `ena_validation_attempts` on `full_seqid` (the assembly prefix plus the annotation
version, e.g. `OG82.ilmn.240313.getorg1770.emma102`) and gives the annotation its own column after
`code`. It rebuilds the table, so rows written before it carry `full_seqid = assembly_prefix` and a
NULL `annotation` until those assemblies are validated again.

**015** adds `ena_submissions`, the submission ledger, and the `ena_submission_status` view. Nothing
in this pipeline writes or reads either one: validation must not depend on submission state. The
table is written by the downstream submitter, keyed on `(full_seqid, webin_mode)`; see
[ena_package.md](../docs/ena_package.md#recording-submissions). `ena_submission_status` joins the
latest validation attempt per `full_seqid` to that ledger.

**016** restores `mitogenome_data.og_num` as a stored generated column,
`(SUBSTRING(og_id FROM 3))::integer`, matching `sample`, `draft_genomes`, `lca` and
`lca_raw_results`. PostgreSQL cannot convert a column in place, so it drops and re-adds it. Never
name `og_num` in an `INSERT` or `UPDATE` column list; a generated column rejects it. The same holds
for the `og_num` that 020 adds to `ena_validation_attempts`.

**017** records the content-addressing of `lca` and `lca_raw_results`: the `content_hash` and
`taxon_rank_db` columns, the `lca_set_content_hash()` trigger function and its two triggers, the
swap of both unique keys to `content_hash`, and the rename of both foreign keys off the ambiguous
`fk_mitogenome`. It is a no-op on the live database, where the change was applied by hand, and is
what lets a rebuilt database satisfy the `ON CONFLICT ON CONSTRAINT lca_content_unique` /
`lca_raw_results_content_unique` targets in `push_lca_blast_results.py` and
`push_lca_raw_results.py`. `content_hash` is trigger-maintained: never name it in an `INSERT` or
`UPDATE` column list.

**018** puts `og_num` back as column 1 of `mitogenome_data`. PostgreSQL cannot reorder columns, so
it rebuilds the table, copies the rows, re-points the three inbound foreign keys and re-grants
`readonly`. The `_1` suffixes on `mitogenome_data_pkey_1`, `mitogenome_data_unique_1` and
`mitogenome_data_depth_method_idx_1` are deliberate: the unsuffixed names belong to the
`mitogenome_data_SS260818` snapshot. **Run it only when no pipeline is writing to
`mitogenome_data`**: it is one transaction, but the foreign-key swap will block or fail against an
in-flight upload. 020 is the same kind of rebuild for `ena_validation_attempts`.
