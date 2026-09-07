/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data upload modules
include { PUSH_MTDNA_ASSM_RESULTS       } from '../../../modules/local/upload_results/mtdna'
include { SPECIES_VALIDATION            } from '../../../modules/local/species_validation'
include { PUSH_MTDNA_ANNOTATION_RESULTS } from '../../../modules/local/upload_results/emma'
include { PUSH_LCA_BLAST_RESULTS        } from '../../../modules/local/upload_results/lca'
include { PUSH_LCA_RAW_RESULTS          } from '../../../modules/local/upload_results/lca_raw'
include { PUSH_ENA_VALIDATION_RESULTS   } from '../../../modules/local/upload_results/ena_validation'
include { PUSH_QC_VALIDATOR             } from '../../../modules/local/upload_results/qc_validator'
include { UPLOAD_RESULTS_SUMMARY        } from '../../../modules/local/upload_results/summary'
include { EVALUATE_QC_CONDITIONS        } from '../../../modules/local/evaluate_qc_conditions'
include { QC_SUMMARY                    } from '../../../modules/local/multiqc/qc_summary'

// Helper functions
include { softwareVersionsToYAML        } from '../../nf-core/utils_nfcore_pipeline'

// Group a per-region result channel ([meta, file], up to one item per annotated
// region) into one [meta, [files]] item per sample, WITHOUT waiting for the
// upstream channel to close.
//
// INVARIANT for this pipeline: a groupTuple whose source is task-derived MUST carry a
// groupKey (or an explicit size:). A bare one is a whole-run barrier by definition -- it
// cannot emit anything until its source channel is complete, which for the LCA/BLAST
// results means every such task in the whole run has finished. That is what stopped a
// finished sample from reaching species validation and QC until the slowest sample in the
// run caught up. groupKey attaches the expected item count to each key, so groupTuple can
// close each sample's group as soon as that sample's own regions are in.
//
// The same operator has since bitten twice more, both times as "one queued task and nothing
// downstream runs at all": the assembly-upload selection in workflows/oceangenomesmitogenomes.nf
// (fixed by disambiguating upstream so no grouping is needed) and ch_ena_validation_inputs in
// subworkflows/local/mitogenome_qc (left in place deliberately -- see the note there for why a
// constant size would be wrong). Where a correct count is not available, remove the need to
// group rather than guessing a size.
//
// The count comes from region_counts ([meta, 0..3]), emitted once per sample by
// MITOGENOME_ANNOTATION. It is exact for both channels because BLAST_BLASTN emits
// one filtered result per region unconditionally, and LCA now emits one lca and
// one lca_raw per region unconditionally too (header-only when a region had no
// valid hits) -- see modules/local/LCA.
//
// combine(by: 0) rather than join(by: 0): join consumes one item per key, so it
// would match only the first of a sample's up to three regions.
// Attach each sample's circularity-check evidence to its QC-gate input, keyed on
// mt_assembly_prefix, WITHOUT waiting for the evidence channel to close.
//
// qc_pairs is [ meta, blast_filtered, annotation_stats ]; circularity_evidence is
// [ mt_assembly_prefix, circularity_check.tsv ], one row per canonical assembly.
//
// A plain join emits each key the moment both sides hold it, so a sample crosses the gate as
// soon as its OWN work is done. The two alternatives both reintroduce a whole-run wait:
// collecting the evidence (toList/collect) yields a value channel that emits only when its
// source CLOSES, and join(..., remainder: true) can only classify an unmatched row at close
// too. Either way one unfinished assembly anywhere in the run pins every finished sample.
//
// The key is mt_assembly_prefix, not the whole meta map: meta gains `circular` between the
// assembly stage and here (and gains it again on the collapse path), and on a -resume the two
// sides carry meta restored from different cache entries, so whole-map equality is not a
// reliable join key -- the same failure documented for the oatk reference join.
//
// Both sides now key on the assembly's IDENTITY, the curated FASTA basename. They did not
// always: the evidence side was keyed by the sample-level assembly prefix while this side had
// been re-stamped to the basename, so for every curated assembly (reseed / _rgj / _collapsed)
// the join simply never matched. A plain join has no way to report that -- it emits nothing
// and says nothing -- so 23 of 168 finished assemblies vanished between species validation and
// the QC gate with no error, no warning and no empty output to notice. Hence the tee below.
def attachCircularityEvidence(qc_pairs, circularity_evidence) {
    def keyed = qc_pairs.map { meta, blast, stats -> [ meta.mt_assembly_prefix, meta, blast, stats ] }

    // Diagnostic tee ONLY. The main path keeps its plain join, because emitting each sample the
    // moment its own work lands is load-bearing here (see above). This branch is allowed the
    // remainder join's whole-run wait precisely because nothing depends on it: it exists to
    // name, at end of run, any assembly that reached the gate and found no evidence to pair
    // with. A silent key miss is what made the original defect invisible; this makes the next
    // one say so.
    keyed
        .map { prefix, _meta, _blast, _stats -> [ prefix, true ] }
        .join(circularity_evidence.map { prefix, _ev -> [ prefix, true ] }, by: 0, remainder: true)
        .filter { items -> items[1] != null && (items.size() < 3 || items[2] == null) }
        .view { items ->
            "WARNING: assembly '${items[0]}' reached the QC gate with no circularity evidence " +
            "under that name and was dropped. Its evidence is keyed by a different name, which " +
            "means an assembly identity was not stamped from its FASTA basename."
        }

    return keyed
        .join(circularity_evidence, by: 0)
        .map { _prefix, meta, blast, stats, evidence -> [ meta, blast, stats, evidence ] }
}

// Gate a channel on each sample's committed assembly upload receipt, keyed on
// mt_assembly_prefix. Returns the rows unchanged, but only once their mitogenome_data row
// exists.
//
// upload_rows is the raw PUSH_MTDNA_ASSM_RESULTS.out.upload, i.e. [ meta, receipt ].
//
// THREE call sites, and they are not all about depth. mitogenome_data carries three inbound
// foreign keys (sql/018_mitogenome_data_og_num_first.sql:215-223): fk_mitogenome_lca on lca,
// fk_mitogenome_lca_raw_results on lca_raw_results, and fk_mitogenome_lca_validation on
// lca_validation. Each is (og_id, tech, seq_date, code) -- exactly mt_assembly_prefix, which
// is why the 4-part prefix is not a compromise key here but the constraint itself. None of
// SPECIES_VALIDATION, PUSH_LCA_RAW_RESULTS or PUSH_LCA_BLAST_RESULTS had a data dependency on
// the parent push, so Nextflow was free to schedule them first, and did: the child insert
// then failed on the foreign key while the run reported success. PUSH_LCA_BLAST_RESULTS is
// covered transitively because it consumes SPECIES_VALIDATION.out.full, so gating
// grouped_blast_lca gates it too -- no fourth call site.
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
// MITOGENOME_QC still ran zero times. The prefix is 1:1 on both sides, so no fan-out.
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

def groupResultsByRegionCount(results, region_counts) {
    return results
        .combine(region_counts, by: 0)
        .map { meta, result_file, n_regions -> [ groupKey(meta, n_regions), result_file ] }
        .groupTuple()
        // Unwrap the GroupKey back to the plain meta map so key equality still
        // holds for the joins downstream.
        .map { key, files -> [ key.getGroupTarget(), files.flatten() ] }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DATE UPLOAD AND SPECIES CHECK WORKFLOW

        - Specific OceanOmics code using the PostgreSQL database.
        - Check the LCA results against nominal species ID and push results to SQL db
        - Determines species that have validated species ID to proceed with QC to prepare
          the sample for submittion to Genbank.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UPLOAD_RESULTS {

    take:
    assembly_results // tuple val(meta), path(fasta), path(assembly_log), path(mito_depth.tsv)
    annotation_results
    blast_filtered_results
    lca_results
    lca_raw_results
    circularity_evidence // tuple val(mt_assembly_prefix), path(circularity_check.tsv) — one row per canonical assembly, both assemblers
    region_counts // tuple val(meta), val(0..3) — CO1/12s/16s regions this sample annotated
    sql_config // params.sql_config

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 
    //
    // MODULE: Pulling out the statistics from the assembly files and updating the SQL database
    //

    PUSH_MTDNA_ASSM_RESULTS (
        assembly_results,
        params.sql_config
    )

    //
    // Re grouping the CO1, 12s and 16s BLAST and LCA results per sample
    //

    grouped_lca   = groupResultsByRegionCount(lca_results, region_counts)
    grouped_blast = groupResultsByRegionCount(blast_filtered_results, region_counts)

    // A sample whose annotation yielded none of CO1/12s/16s never enters BLAST or
    // LCA, so it appears in neither grouped channel. Without this it would be
    // dropped by the join below and never reach SPECIES_VALIDATION at all -- no
    // QC verdict, no lca_results row, silently absent from the run's output.
    // Feed it through on empty stand-ins instead: species_validation.py combines
    // them to a header-only lca_combined and an empty blast_combined, finds no
    // nominal-species hit, and records the sample as not validated. That is the
    // correct outcome for a sample with nothing to validate against, and it is
    // reached immediately rather than at the end of the run.
    def empty_lca_file   = file("${projectDir}/assets/placeholders/empty_lca.tsv", checkIfExists: true)
    def empty_blast_file = file("${projectDir}/assets/placeholders/empty_blast_filtered.tsv", checkIfExists: true)

    ch_zero_region_blast_lca = region_counts
        .filter { _meta, n_regions -> n_regions == 0 }
        .map { meta, _n_regions -> [ meta, [ empty_blast_file ], [ empty_lca_file ] ] }

    grouped_blast_lca = grouped_blast
        .join(grouped_lca, by: 0)
        .mix(ch_zero_region_blast_lca)

    // Gate on the committed mitogenome_data row before species validation writes
    // lca_validation (and, transitively, before PUSH_LCA_BLAST_RESULTS writes lca):
    // both carry a foreign key to it. Note the zero-region stand-ins are mixed in ABOVE
    // this line on purpose -- they arrive from region_counts rather than from the grouped
    // channels, and they still get an lca_validation row, so they need gating too.
    grouped_blast_lca_gated = gateOnAssemblyReceipt(
        grouped_blast_lca,
        PUSH_MTDNA_ASSM_RESULTS.out.upload,
        'species validation'
    )

    //
    // MODULE: Checking the LCA results against the nominal species ID in the SQL database
    //

    SPECIES_VALIDATION (
        grouped_blast_lca_gated, // tuple val(meta), path(blast_filtered), path(lca_filtered)
        sql_config // params.sql_config
    )

    //
    // MODULE: Calculating the statistics of the annotations and updating the SQL database
    //

    PUSH_MTDNA_ANNOTATION_RESULTS (
        annotation_results, // tuple val(meta), path("annotation/*")
        sql_config // params.sql_config
    )

    //
    // MODULE: Updating the SQL database with the LCA and filtered BLAST results
    //

    PUSH_LCA_BLAST_RESULTS (
        SPECIES_VALIDATION.out.full, // tuple path ("lca_combined.${mt_assembly_prefix}.tsv"), path ("blast_combined.${mt_assembly_prefix}.tsv"),
        sql_config // params.sql_config
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
    // only SPECIES_VALIDATION would close one FK child and leave this one open.
    grouped_lca_raw_gated = gateOnAssemblyReceipt(
        grouped_lca_raw,
        PUSH_MTDNA_ASSM_RESULTS.out.upload,
        'raw LCA upload'
    )

    PUSH_LCA_RAW_RESULTS (
        grouped_lca_raw_gated, // tuple val(meta), path(lca_raw.*.tsv)
        sql_config // params.sql_config
    )

    //
    // SUBWORKFLOW: Evaluate the mitogenome and if it can continue on to final QC
    //
    /*  This solution provides a conditional QC subworkflow that evaluates two conditions and only proceeds with QC processes when **both** are satisfied:
        1. ✅ **LCA results**: Any row has the nominal species ID in the filtered blast
            results: `Found_in_blast_YN = "Yes"`
        2. ✅ **Annotation CSV**: The `passed` column value is `"yes"` when there are no
            missing genes and the genes are in the right order.
    */
    
    // Attach the per-sample circularity-check evidence so the gate can block length/repeat
    // anomalies (concatemers, control-region VNTRs, unresolved over-length assemblies) from
    // progressing to QC. Assemblies that never reached a check (failed / no-contig /
    // precomputed with no *_check.tsv on disk) carry the empty-check stand-in from the
    // assembly stage, so they are never dropped and never blocked on this condition.
    //
    // The evidence arrives as one row per canonical assembly, keyed by mt_assembly_prefix,
    // so it can be attached with a plain join: each sample crosses the gate as soon as its
    // OWN evidence, species validation and annotation stats are in.
    //
    // This previously collected the whole evidence channel into a lookup map
    // (circularity_evidence.toList()) to avoid a remainder join's whole-run wait. That traded
    // one barrier for another: toList() is a value channel that emits only when its source
    // CLOSES, and the source runs back to the assembly subworkflows, so a single assembly
    // still queued held every finished sample at this gate. Run
    // mitogenomes-missing-audit-5 hit exactly that -- 135 samples validated and uploaded,
    // EVALUATE_QC_CONDITIONS never ran once, because one sample's REFERENCE_RANK was stuck
    // behind a maintenance reservation.
    //
    // The fix is upstream totality rather than a cleverer operator here: every emitted
    // assembly carries an evidence row (assets/placeholders/empty_circularity_check.tsv where no check
    // ran), so there is no "missing evidence" case left for a lookup fallback or a remainder
    // path to cover, and no reason to wait for the channel to close. Deliberately no
    // placeholder fallback at this point: if the contract upstream ever breaks, the sample
    // should visibly fail to reach QC rather than be silently gated on a stand-in that says
    // "no anomaly".
    ch_qc_conditions = attachCircularityEvidence(
        SPECIES_VALIDATION.out.summary.join(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats, by: 0),
        circularity_evidence
    )

    //
    // MODULE: evaluating the results to determine if to process the sample through QC
    //

    EVALUATE_QC_CONDITIONS (
        ch_qc_conditions // tuple val(meta), path(blast_filtered), path(annotation_stats.csv), path(circularity_check.tsv)
    )

    // Filter for samples that meet both conditions, then gate each on its own committed
    // assembly upload row. See gateOnAssemblyReceipt for why the key is the prefix:
    // when this was a whole-meta join it matched NOTHING once the assembly-upload barrier was
    // removed -- 123 samples cleared the gate, all 123 had a matching upload row by prefix,
    // and MITOGENOME_QC still received zero.
    ch_qc_ready = gateOnAssemblyReceipt(
        EVALUATE_QC_CONDITIONS.out.evaluation
            .map { meta, species_file, proceed_file, circular_file ->
                def species_name = species_file.text.trim()
                def proceed_qc = proceed_file.text.trim()
                def circular = circular_file.text.trim()
                return [ meta, species_name, proceed_qc, circular ]
            }
            .filter { meta, species_name, proceed_qc, circular ->
                proceed_qc == "true"
            },
        PUSH_MTDNA_ASSM_RESULTS.out.upload,
        'the QC gate'
    )

    // Log samples that will proceed to QC
    ch_qc_ready.view { meta, species_name, proceed_qc, circular ->
        "Sample ${meta.id} will proceed to QC with species: ${species_name} (circular: ${circular})"
    }

    // Filter for samples that dont meet the conditions
    ch_not_qc_ready = EVALUATE_QC_CONDITIONS.out.evaluation
        .join(EVALUATE_QC_CONDITIONS.out.reason, by: 0)
        .map { meta, species_file, proceed_file, circular_file, reason_file ->
            def species_name = species_file.text.trim()
            def proceed_qc = proceed_file.text.trim()
            def circular = circular_file.text.trim()
            def held_reason = reason_file.text.trim()
            return [ meta, species_name, proceed_qc, circular, held_reason ]
        }
        .filter { meta, species_name, proceed_qc, circular, held_reason ->
            proceed_qc == "false"
        }
        .view { meta, species_name, proceed_qc, circular, held_reason ->
            "Sample ${meta.id} will NOT proceed to QC - ${held_reason ?: 'conditions not met'}"
        }

    // Headerless per-sample fragments for the run-level held_samples.tsv. These
    // samples are filtered out before MITOGENOME_QC, so this subworkflow is the
    // only place their hold is recorded.
    ch_held_fragments = ch_not_qc_ready
        .collectFile { meta, _species, _proceed, _circular, held_reason ->
            [ "${meta.mt_assembly_prefix}.held.tsv",
              "${meta.id}\t${meta.mt_assembly_prefix}\tPRE_QC\tproceed_qc=false: ${held_reason ?: 'conditions not met'}\n" ]
        }


    //
    // Build a per-sample QC summary TSV for MultiQC
    qc_summary_input = EVALUATE_QC_CONDITIONS.out.evaluation
        .join(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats, by: 0)
        .map { meta, species_file, proceed_file, circular_file, annotation_csv -> [ meta, species_file, proceed_file, annotation_csv ] }
    QC_SUMMARY (
        qc_summary_input // tuple val(meta), path(species_name.txt), path(proceed_qc.txt), path(annotation_stats.csv)
    )

    
    //
    // MODULE: Consolidate the per-step SQL upload status files into a single
    //         summary TSV + detailed appendix so all upload outcomes can be
    //         reviewed in one place instead of grep'ing many small files.
    //

    ch_upload_status_files = Channel.empty()
        .mix(PUSH_MTDNA_ASSM_RESULTS.out.upload.map { _meta, file -> file })
        .mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.upload.map { _meta, f -> f })
        .mix(SPECIES_VALIDATION.out.upload.map { _meta, f -> f })
        .mix(PUSH_LCA_BLAST_RESULTS.out.upload)
        .mix(PUSH_LCA_RAW_RESULTS.out.upload)

    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    //  - Species validation outputs (per-sample TSVs)
    //  - Annotation statistics CSVs
    //  - QC evaluation flags (for quick visibility in report) and summary table
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_VALIDATION.out.summary.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.evaluation.map { meta, species_file, proceed_file, circular_file -> proceed_file })
    ch_multiqc_files = ch_multiqc_files.mix(QC_SUMMARY.out.table.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ASSM_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_VALIDATION.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_LCA_BLAST_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_LCA_RAW_RESULTS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(PUSH_MTDNA_ASSM_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(SPECIES_VALIDATION.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_BLAST_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_RAW_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(EVALUATE_QC_CONDITIONS.out.versions)



    //
    // Emit outputs
    //

    emit:
    qc_ready    = ch_qc_ready                   // channel: [ val(meta), val(species_name), val(proceed_qc true/false), val(circular true/false) ]
    held_fragments = ch_held_fragments         // channel: path(<prefix>.held.tsv) — one row per pre-QC hold
    assembly_summary_files = PUSH_MTDNA_ANNOTATION_RESULTS.out.stats.map { meta, stats -> stats }
    upload_status_files = ch_upload_status_files
    multiqc_files = ch_multiqc_files            // channel: [ path(multiqc_files) ]
    versions = ch_versions             // channel: [ path(versions.yml) ]
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
