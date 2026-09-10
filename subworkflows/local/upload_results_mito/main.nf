/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS

        Every process in this file is a PUSHER: it writes to the OceanOmics
        PostgreSQL database and does nothing else. Anything that computes a verdict
        or a report now lives in subworkflows/local/mitogenome_qc, which runs
        unconditionally. See the note at the top of that file for why the two were
        separated.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data upload modules
include { PUSH_MTDNA_ASSM_RESULTS       } from '../../../modules/local/upload_results/mtdna'
include { PUSH_SPECIES_VALIDATION       } from '../../../modules/local/upload_results/species_validation'
include { PUSH_MTDNA_ANNOTATION_RESULTS } from '../../../modules/local/upload_results/emma'
include { PUSH_LCA_BLAST_RESULTS        } from '../../../modules/local/upload_results/lca'
include { PUSH_LCA_RAW_RESULTS          } from '../../../modules/local/upload_results/lca_raw'
include { PUSH_ENA_VALIDATION_RESULTS   } from '../../../modules/local/upload_results/ena_validation'
include { PUSH_QC_VALIDATOR             } from '../../../modules/local/upload_results/qc_validator'
include { UPLOAD_RESULTS_SUMMARY        } from '../../../modules/local/upload_results/summary'

// Helper functions
include { softwareVersionsToYAML        } from '../../nf-core/utils_nfcore_pipeline'
// Grouping helper shared with the QC subworkflow that owns it. The dependency runs
// this way round on purpose: uploads depend on QC, never the reverse.
include { groupResultsByRegionCount     } from '../mitogenome_qc/main'

// Gate a channel on each sample's committed assembly upload receipt, keyed on
// mt_assembly_prefix. Returns the rows unchanged, but only once their mitogenome_data row
// exists.
//
// upload_rows is the raw PUSH_MTDNA_ASSM_RESULTS.out.upload, i.e. [ meta, receipt ].
//
// ONE CALL SITE PER PUSHER, and they are not all about depth. mitogenome_data carries three
// inbound foreign keys (sql/018_mitogenome_data_og_num_first.sql:215-223): fk_mitogenome_lca
// on lca, fk_mitogenome_lca_raw_results on lca_raw_results, and fk_mitogenome_lca_validation
// on lca_validation. Each is (og_id, tech, seq_date, code) -- exactly mt_assembly_prefix,
// which is why the 4-part prefix is not a compromise key here but the constraint itself. None
// of these pushers has a data dependency on the parent push, so Nextflow is free to schedule
// them first, and did: the child insert then failed on the foreign key while the run reported
// success.
//
// Every pusher is now gated EXPLICITLY at its own call site. Two of them used to inherit the
// gate transitively instead, through a SPECIES_VALIDATION whose input was gated:
// PUSH_LCA_BLAST_RESULTS consumed SPECIES_VALIDATION.out.full, and the annotation push
// consumed a channel joined against it. That worked only while SPECIES_VALIDATION was itself
// a writer sitting inside this subworkflow. It is now DB-free QC that runs unconditionally
// upstream (subworkflows/local/mitogenome_qc), so nothing it emits is ordered against the
// parent row any more, and an inherited gate would silently be no gate at all. If you add a
// pusher, gate it here; do not assume its inputs were gated for you.
//
// The `annotation` component of the 5-part lca_validation identity is deliberately NOT part
// of this key: it is not in meta at this stage (species_validation.py derives it inside the
// task from the BLAST query id), and the FK does not use it. Anything that identifies or
// repairs an lca_validation ROW still needs all five parts.
//
// The receipt is a pure ordering dependency: ENA metadata reads mean_depth from
// mitogenome_data, so QC must not race ahead of the committed row. The key must be the prefix
// and not the whole meta map -- exactly the lesson stated for attachCircularityEvidence above,
// which this join was ignoring. The two sides arrive down different lineages and on a -resume
// are restored from separate cache entries, so a whole-map join matches nothing and drops
// every sample SILENTLY instead of failing -- which is what it did the moment the
// assembly-upload barrier upstream was removed and this became the binding constraint:
// 123 samples cleared the gate, all 123 had a matching upload row by prefix, and
// ENA_SUBMISSION_PREP still ran zero times. The prefix is 1:1 on both sides, so no fan-out.
//
// Both sides key on the assembly's IDENTITY, and the upload side is now built from the
// SANITISED assembly, so the receipt exists under the same curated name the QC row carries.
// It did not used to: the canonical upload row was filed under the assembly-stage name, so a
// _collapsed or _concat assembly had no receipt under its own name and would have been dropped
// here even after the evidence join was fixed. reseed / _rgj happened to survive only because
// they get a provenance row of their own that collided with the right name by accident.
def gateOnAssemblyReceipt(rows, upload_rows, label) {
    // rows is any tuple whose FIRST element is meta; the shape is preserved. Keying by
    // index rather than destructuring is what lets one helper serve tuples of different
    // arity -- [meta, blast, lca], [meta, lca_raw] and [meta, species, proceed, circular]
    // all gate the same way, and hand-rolling a third copy is how the REF_GENES drift
    // documented in bin/mito_gene_order.py happened.
    def keyed = rows.map { row -> [ (row[0].mt_assembly_prefix) ] + row }
    def receipts = upload_rows.map { meta, receipt -> [ meta.mt_assembly_prefix, receipt ] }

    // Diagnostic tee ONLY, same shape and same reasoning as attachCircularityEvidence:
    // a plain join that misses its key emits nothing and says nothing, and introducing a
    // new join is introducing a new silent-drop risk. This branch is allowed the remainder
    // join's whole-run wait precisely because nothing depends on it.
    keyed
        .map { items -> [ items[0], true ] }
        .join(receipts.map { prefix, _receipt -> [ prefix, true ] }, by: 0, remainder: true)
        .filter { items -> items[1] != null && (items.size() < 3 || items[2] == null) }
        .view { items ->
            "WARNING: assembly '${items[0]}' reached ${label} with no committed " +
            "mitogenome_data upload receipt under that name. Its FK parent row was never " +
            "written, or the receipt is keyed by a different name."
        }

    return keyed
        .join(receipts, by: 0)
        // Drop the join key we added at the front and the receipt the join appended at the
        // back, restoring the caller's original tuple shape.
        .map { items -> items[1..-2] }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DATE UPLOAD AND SPECIES CHECK WORKFLOW

        - Specific OceanOmics code using the PostgreSQL database.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UPLOAD_RESULTS {

    take:
    assembly_results        // tuple val(meta), path(fasta), path(assembly_log), path(mito_depth.tsv)
    annotation_stats        // tuple val(meta), path(annotation_stats.csv) — from MITOGENOME_QC
    validation_records      // tuple val(meta), path(validation_record.json) — from MITOGENOME_QC
    species_validation_full // tuple val(meta), path(lca_combined), path(blast_combined) — from MITOGENOME_QC
    lca_raw_results
    region_counts           // tuple val(meta), val(0..3) — CO1/12s/16s regions this sample annotated
    sql_config              // params.sql_config

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Pulling out the statistics from the assembly files and updating the SQL database
    //
    // This is the FK parent for everything below, so it is first and everything else
    // is gated on its receipt.
    //

    PUSH_MTDNA_ASSM_RESULTS (
        assembly_results,
        params.sql_config
    )

    //
    // MODULE: Writing the lca_validation row from the record MITOGENOME_QC produced.
    //         lca_validation carries fk_mitogenome_lca_validation, so it is gated.
    //

    PUSH_SPECIES_VALIDATION (
        gateOnAssemblyReceipt(
            validation_records,
            PUSH_MTDNA_ASSM_RESULTS.out.upload,
            'species validation upload'
        ),
        sql_config
    )

    //
    // MODULE: Updating the SQL database with the annotation statistics
    //
    // Gated explicitly. It upserts into mitogenome_data itself rather than a child
    // table, so this is an ordering constraint rather than a foreign key: without it
    // the annotation columns can be written first and the assembly push can then
    // overwrite the row it finds. It used to inherit the ordering by consuming a
    // channel joined against the gated SPECIES_VALIDATION output; that inheritance is
    // gone now that species validation is DB-free QC upstream.
    //

    PUSH_MTDNA_ANNOTATION_RESULTS (
        gateOnAssemblyReceipt(
            annotation_stats,
            PUSH_MTDNA_ASSM_RESULTS.out.upload,
            'annotation upload'
        ),
        sql_config
    )

    //
    // MODULE: Updating the SQL database with the LCA and filtered BLAST results
    //
    // lca carries fk_mitogenome_lca. Gated explicitly for the same reason as the
    // annotation push: it used to be covered transitively through SPECIES_VALIDATION.
    //

    PUSH_LCA_BLAST_RESULTS (
        gateOnAssemblyReceipt(
            species_validation_full, // tuple val(meta), path(lca_combined), path(blast_combined)
            PUSH_MTDNA_ASSM_RESULTS.out.upload,
            'LCA/BLAST upload'
        ),
        sql_config
    )

    //
    // MODULE: Pushing raw per-hit LCA rows into the lca_raw_results table.
    //         Group the per-region lca_raw files by sample so each sample gets
    //         one upload call with all its regions.
    //
    //         Zero-region samples are deliberately absent here: they have no raw
    //         hits to insert, so there is no row to push.
    grouped_lca_raw = groupResultsByRegionCount(lca_raw_results, region_counts)

    // lca_raw_results carries fk_mitogenome_lca_raw_results, so the same receipt gate
    // applies. This is the LARGER of the two failure populations, not the smaller: gating
    // only the validation upload would close one FK child and leave this one open.
    grouped_lca_raw_gated = gateOnAssemblyReceipt(
        grouped_lca_raw,
        PUSH_MTDNA_ASSM_RESULTS.out.upload,
        'raw LCA upload'
    )

    PUSH_LCA_RAW_RESULTS (
        grouped_lca_raw_gated, // tuple val(meta), path(lca_raw.*.tsv)
        sql_config
    )

    //
    // MODULE: Consolidate the per-step SQL upload status files into a single
    //         summary TSV + detailed appendix so all upload outcomes can be
    //         reviewed in one place instead of grep'ing many small files.
    //

    ch_upload_status_files = Channel.empty()
        .mix(PUSH_MTDNA_ASSM_RESULTS.out.upload.map { _meta, file -> file })
        .mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.upload.map { _meta, f -> f })
        .mix(PUSH_SPECIES_VALIDATION.out.upload.map { _meta, f -> f })
        .mix(PUSH_LCA_BLAST_RESULTS.out.upload)
        .mix(PUSH_LCA_RAW_RESULTS.out.upload)

    //
    // Subworkflow finishing steps.
    //

    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ASSM_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_SPECIES_VALIDATION.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_LCA_BLAST_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_LCA_RAW_RESULTS.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(PUSH_MTDNA_ASSM_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_SPECIES_VALIDATION.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_BLAST_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_RAW_RESULTS.out.versions.first())

    //
    // Emit outputs
    //

    emit:
    // The committed mitogenome_data receipts, the FK parent for every other write in
    // here. No consumer OUTSIDE this subworkflow uses them any more: they were exposed
    // to order ENA_SUBMISSION_PREP behind the committed row, back when ENA metadata read
    // mean_depth out of it. Submission prep now gets that number from the run's own depth
    // TSV and needs no ordering at all. Kept exported because the receipt stream is what
    // tests/upload_receipt_join exercises, and because a future pusher outside this
    // subworkflow would need exactly this to gate on.
    assembly_receipts   = PUSH_MTDNA_ASSM_RESULTS.out.upload
    upload_status_files = ch_upload_status_files
    multiqc_files       = ch_multiqc_files // channel: [ path(multiqc_files) ]
    versions            = ch_versions      // channel: [ path(versions.yml) ]
}
/*
 * Upload the post-QC ENA validation records and compile the final SQL upload
 * report only after all biological and submission-readiness gates have run.
 */
workflow UPLOAD_ENA_RESULTS {

    take:
    ena_validation_records // tuple val(meta), path(*.ena_validation_result.tsv)
    prior_upload_status_files
    sql_config

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    PUSH_ENA_VALIDATION_RESULTS(ena_validation_records, sql_config)

    // Second species-ID validator. A sample that cleared every QC gate
    // (submission_ready = true) has been checked harder by the pipeline than a
    // second reviewer would manage, so the pipeline signs off as
    // lca_validation.validator_2. Independent of the call above -- different
    // table, no ordering dependency -- and it never overwrites a validator_2
    // that a human already filled in.
    PUSH_QC_VALIDATOR(ena_validation_records, sql_config)

    ch_all_upload_status_files = prior_upload_status_files
        .mix(PUSH_ENA_VALIDATION_RESULTS.out.upload.map { _meta, upload -> upload })
        .mix(PUSH_QC_VALIDATOR.out.upload.map { _meta, upload -> upload })

    UPLOAD_RESULTS_SUMMARY(ch_all_upload_status_files.collect())

    ch_multiqc_files = ch_multiqc_files
        .mix(PUSH_ENA_VALIDATION_RESULTS.out.tool_params.collect { it[1] })
        .mix(PUSH_QC_VALIDATOR.out.tool_params.collect { it[1] })
        .mix(UPLOAD_RESULTS_SUMMARY.out.summary)
        .mix(UPLOAD_RESULTS_SUMMARY.out.tool_params)
    ch_versions = ch_versions
        .mix(PUSH_ENA_VALIDATION_RESULTS.out.versions.first())
        .mix(PUSH_QC_VALIDATOR.out.versions.first())
        .mix(UPLOAD_RESULTS_SUMMARY.out.versions)

    emit:
    upload_logs = PUSH_ENA_VALIDATION_RESULTS.out.upload
    upload_summary = UPLOAD_RESULTS_SUMMARY.out.summary
    upload_appendix = UPLOAD_RESULTS_SUMMARY.out.appendix
    multiqc_files = ch_multiqc_files
    versions = ch_versions
}
