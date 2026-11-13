/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MITOGENOME QC

        - Decides whether an assembly is good enough to proceed: the species
          comparison, the annotation statistics, the gate that combines them, and
          the per-sample QC summary.

        - Contains NO database access, on purpose. Every one of these steps used to
          live inside UPLOAD_RESULTS behind
          `if (!params.skip_upload_results && params.sql_config)`, which meant the
          only way to learn whether a mitogenome was any good was to also write it to
          PostgreSQL. A run with --skip_upload_results, or with no --sql_config at
          all -- which is how the invertebrate work is developed -- produced no gene
          counts, no completeness verdict, no held-samples accounting, and an
          assembly summary whose annotation columns were empty for every sample.

        - The rule this file exists to enforce: a process either produces a verdict
          (here) or writes to the database (upload_results_mito), never both.
          SPECIES_VALIDATION and the old PUSH_MTDNA_ANNOTATION_RESULTS each did both,
          and that is the whole reason QC was reachable only through the uploader.

        - Runs downstream of annotation and upstream of ENA_SUBMISSION_PREP, which
          takes the samples this gate releases and prepares them for submission.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPECIES_VALIDATION     } from '../../../modules/local/species_validation'
include { ANNOTATION_STATS       } from '../../../modules/local/annotation_stats'
include { EVALUATE_QC_CONDITIONS } from '../../../modules/local/evaluate_qc_conditions'
include { QC_SUMMARY             } from '../../../modules/local/multiqc/qc_summary'

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
// subworkflows/local/ena_submission_prep (left in place deliberately -- see the note there for why a
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
def groupResultsByRegionCount(results, region_counts) {
    return results
        .combine(region_counts, by: 0)
        .map { meta, result_file, n_regions -> [ groupKey(meta, n_regions), result_file ] }
        .groupTuple()
        // Unwrap the GroupKey back to the plain meta map so key equality still
        // holds for the joins downstream.
        .map { key, files -> [ key.getGroupTarget(), files.flatten() ] }
}

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



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME QC WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_QC {

    take:
    annotation_results     // tuple val(meta), path("annotation/*")
    blast_filtered_results // tuple val(meta), path(blast filtered) — one per annotated region
    lca_results            // tuple val(meta), path(lca) — one per annotated region
    circularity_evidence   // tuple val(mt_assembly_prefix), path(circularity_check.tsv) — one row per canonical assembly
    region_counts          // tuple val(meta), val(0..3) — CO1/12s/16s regions this sample annotated

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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

    //
    // MODULE: Checking the LCA results against the nominal species ID
    //
    // No upload-receipt gate here, and none is needed: this writes nothing to the
    // database. The gate exists to order FK children after their mitogenome_data
    // parent, and the only FK child in this lineage is now PUSH_SPECIES_VALIDATION,
    // which is gated at its own call site in upload_results_mito.
    //
    SPECIES_VALIDATION (
        grouped_blast_lca // tuple val(meta), path(blast_filtered), path(lca_filtered)
    )

    //
    // MODULE: Calculating the statistics of the annotations
    //

    // Attach each assembly's lca_combined so annotation_stats.py can record whether
    // the BLAST-derived lineage agreed with the rank a gene-order variant rule
    // matched on. Advisory only -- it holds nothing.
    //
    // A PLAIN join, deliberately. remainder: true is a whole-run barrier here (see
    // the note in attachCircularityEvidence above) and has already broken a run. A
    // plain join is only safe because SPECIES_VALIDATION is TOTAL over
    // annotation_results on both paths: the main path derives region_counts
    // directly from the annotation results, the zero-region stand-in above catches
    // the rest, and the precomputed path was made total by construction in
    // workflows/oceangenomesmitogenomes.nf. Before that fix a precomputed assembly
    // with annotation files but no BLAST files reached neither, and this join would
    // have silently dropped it -- it gets its annotation stats computed today, so
    // that would have been a regression.
    //
    // The tee below is the standing insurance: it names, at end of run, any
    // assembly that reached here with no lca_combined to pair with, because a plain
    // join that misses emits nothing and says nothing.
    ch_annotation_lca = SPECIES_VALIDATION.out.full
        .map { meta, lca_combined, _blast_combined ->
            [ meta.mt_assembly_prefix, lca_combined ]
        }

    annotation_results
        .map { meta, _files -> [ meta.mt_assembly_prefix, true ] }
        .join(ch_annotation_lca.map { prefix, _f -> [ prefix, true ] }, by: 0, remainder: true)
        .filter { items -> items[1] != null && (items.size() < 3 || items[2] == null) }
        .view { items ->
            "WARNING: assembly '${items[0]}' reached the annotation statistics with no " +
            "lca_combined under that name and was dropped. SPECIES_VALIDATION is not " +
            "total over annotation_results, which is a totality bug upstream, not a " +
            "reason to loosen this join."
        }

    // tuple val(meta), path("annotation/*"), path(lca_combined)
    ch_annotation_with_lca = annotation_results
        .map { meta, files -> [ meta.mt_assembly_prefix, meta, files ] }
        .join(ch_annotation_lca, by: 0)
        .map { _prefix, meta, files, lca_combined -> [ meta, files, lca_combined ] }

    ANNOTATION_STATS (
        ch_annotation_with_lca
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
        SPECIES_VALIDATION.out.summary.join(ANNOTATION_STATS.out.stats, by: 0),
        circularity_evidence
    )

    //
    // MODULE: evaluating the results to determine if to process the sample through QC
    //

    EVALUATE_QC_CONDITIONS (
        ch_qc_conditions // tuple val(meta), path(blast_filtered), path(annotation_stats.csv), path(circularity_check.tsv)
    )

    // Filter for samples that meet both conditions.
    //
    // This used to be gated on the assembly upload receipt, because ENA_SUBMISSION_PREP
    // read mean_depth back out of mitogenome_data and could not race ahead of the
    // committed row. Gating the VERDICT on a database write is what made QC unreachable
    // without a database, which is the defect this subworkflow fixes. The gate first
    // moved to ENA_SUBMISSION_PREP's own input, and is now gone entirely: coverage
    // travels to submission prep as the depth TSV this run measured, so nothing
    // downstream waits on a write. gateOnAssemblyReceipt still guards the pushers.
    ch_qc_ready = EVALUATE_QC_CONDITIONS.out.evaluation
        .map { meta, species_file, proceed_file, circular_file ->
            def species_name = species_file.text.trim()
            def proceed_qc = proceed_file.text.trim()
            def circular = circular_file.text.trim()
            return [ meta, species_name, proceed_qc, circular ]
        }
        .filter { meta, species_name, proceed_qc, circular ->
            proceed_qc == "true"
        }

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
    // samples are filtered out before ENA_SUBMISSION_PREP, so this subworkflow is the
    // only place their hold is recorded.
    //
    // The filename carries the STAGE as well as the prefix. It did not always, and every
    // held source used the same "<prefix>.held.tsv" name: the fragments are collectFile
    // outputs from different subworkflows that are then mixed and staged flat into
    // fragments/, so two genuine holds on one assembly at two stages collided under one
    // name, and COMPILE_HELD_SAMPLES' sort -u then collapsed whatever survived. An
    // assembly can legitimately be held more than once.
    ch_held_fragments = ch_not_qc_ready
        .collectFile { meta, _species, _proceed, _circular, held_reason ->
            [ "${meta.mt_assembly_prefix}.PRE_QC.held.tsv",
              "${meta.id}\t${meta.mt_assembly_prefix}\tPRE_QC\tproceed_qc=false: ${held_reason ?: 'conditions not met'}\n" ]
        }

    //
    // Build a per-sample QC summary TSV for MultiQC
    qc_summary_input = EVALUATE_QC_CONDITIONS.out.evaluation
        .join(ANNOTATION_STATS.out.stats, by: 0)
        .map { meta, species_file, proceed_file, circular_file, annotation_csv -> [ meta, species_file, proceed_file, annotation_csv ] }
    QC_SUMMARY (
        qc_summary_input // tuple val(meta), path(species_name.txt), path(proceed_qc.txt), path(annotation_stats.csv)
    )

    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    //  - Species validation outputs (per-sample TSVs)
    //  - Annotation statistics CSVs
    //  - QC evaluation flags (for quick visibility in report) and summary table
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_VALIDATION.out.summary.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATION_STATS.out.stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.evaluation.map { meta, species_file, proceed_file, circular_file -> proceed_file })
    ch_multiqc_files = ch_multiqc_files.mix(QC_SUMMARY.out.table.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SPECIES_VALIDATION.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATION_STATS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(EVALUATE_QC_CONDITIONS.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(SPECIES_VALIDATION.out.versions.first())
    ch_versions = ch_versions.mix(ANNOTATION_STATS.out.versions.first())
    ch_versions = ch_versions.mix(EVALUATE_QC_CONDITIONS.out.versions)

    //
    // Emit outputs
    //

    emit:
    qc_ready               = ch_qc_ready                  // channel: [ val(meta), val(species_name), val(proceed_qc true/false), val(circular true/false) ]
    held_fragments         = ch_held_fragments            // channel: path(<prefix>.PRE_QC.held.tsv) — one row per pre-QC hold
    // The annotation statistics CSV, for MITOGENOME_ASSEMBLY_SUMMARY. This is the
    // channel whose absence left num_genes/num_cds/missing_genes empty for every
    // sample in any run without a database.
    assembly_summary_files = ANNOTATION_STATS.out.stats.map { _meta, stats -> stats }
    annotation_stats       = ANNOTATION_STATS.out.stats   // channel: [ val(meta), path(annotation_stats.csv) ]
    validation_records     = SPECIES_VALIDATION.out.validation_record // channel: [ val(meta), path(validation_record.json) ]
    species_validation_full = SPECIES_VALIDATION.out.full // channel: [ val(meta), path(lca_combined), path(blast_combined) ]
    evaluation             = EVALUATE_QC_CONDITIONS.out.evaluation
    multiqc_files          = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions               = ch_versions                  // channel: [ path(versions.yml) ]
}
