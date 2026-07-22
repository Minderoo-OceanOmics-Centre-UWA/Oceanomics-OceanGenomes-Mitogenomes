/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Pipeline subworkflows
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { MULTIQC_PER_SAMPLE     } from '../modules/local/multiqc/per_sample'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oceangenomesmitogenomes_pipeline'

// Mitogenome assembly subworkflows
include { MITOGENOME_ASSEMBLY_GETORG       } from '../subworkflows/local/mitogenome_assembly/getorganelle'
include { MITOGENOME_ASSEMBLY_MITOHIFI     } from '../subworkflows/local/mitogenome_assembly/mitohifi'
include { MITOGENOME_ANNOTATION     } from '../subworkflows/local/mitogenome_annotation_lca'
include { COLLAPSE_CONCATEMER       } from '../modules/local/collapse_concatemer'
include { UPLOAD_RESULTS; UPLOAD_ENA_RESULTS } from '../subworkflows/local/upload_results_mito'
include { MITOGENOME_QC             } from '../subworkflows/local/mitogenome_qc'
include { SANITISE_FASTA           } from '../modules/local/sanitise_fasta/main'
include { MITOGENOME_ASSEMBLY_SUMMARY } from '../modules/local/multiqc/mitogenome_assembly_summary'

// Using local module SANITISE_FASTA (see modules/local/sanitise_fasta)

def evidenceAnomalyType(tsv) {
    try {
        def rows = tsv.text.readLines()
        if (rows.size() < 2) return ''
        def header = rows[0].split('\t', -1)
        def idx = header.findIndexOf { it.trim() == 'anomaly_type' }
        def values = rows[1].split('\t', -1)
        return idx >= 0 && idx < values.size() ? values[idx].trim().toLowerCase() : ''
    } catch (ignored) {
        return ''
    }
}

def fastaSequenceLength(fasta) {
    try {
        return fasta.text.readLines().findAll { !it.startsWith('>') }.sum { it.trim().size() } ?: 0
    } catch (ignored) {
        return 0
    }
}

def evidenceFinalVerdictCircular(tsv) {
    try {
        def rows = tsv.text.readLines()
        if (rows.size() < 2) return null
        def header = rows[0].split('\t', -1)
        def idx = header.findIndexOf { it.trim() == 'final_verdict_circular' }
        def values = rows[1].split('\t', -1)
        if (idx < 0 || idx >= values.size()) return null
        def value = values[idx].trim().toLowerCase()
        return value == 'true' ? true : (value == 'false' ? false : null)
    } catch (ignored) {
        return null
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OCEANGENOMESMITOGENOMES {

    take:
    getorg_input    // tuple for getorganelle - tuple(meta, reads)
    mitohifi_input // tuple for mitohifi - tuple(meta, reads)
    curated_blast_db // params.curated_blast_db
    nt_blast_db // params.nt_blast_db
    mitos_refdb // params.mitos_refdb (MITOS2 RefSeq reference data dir)
    sql_config // params.sql_config
    organelle_type // params.organelle_type "animal_mt"

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_assembly_summary_files = Channel.empty()

    // Map samplesheet meta by sample id + sequencing type + date for reuse with precomputed files
    ch_samplesheet_meta = getorg_input
        .map { meta, _reads -> [ [meta.id, meta.sequencing_type, meta.date], meta ] }
        .mix(mitohifi_input.map { meta, _reads -> [ [meta.id, meta.sequencing_type, meta.date], meta ] })
        .groupTuple()
        .map { sample_key, metas -> [ sample_key, metas[0] ] }

    /* The if and else statements in this workflow are for when steps are skipped in the nextflow_run script.
        What it is doing is 'if' this processes isnt skipped then run the subworkflow and provide the standard outputs.
        
        Then 'else if' is if the process has been skipped then we will create the required input for the following
        subworkflows using the predefined file paths for the output in the nextflow.config, the "precomputed_*" file
        paths. 
            The relevant information is extracted to create the meta map from the parts of the file name. This is
            assuming that files are named with the sample id, then the type of sequencing, the date and then the other
            information in the file name, all seperated by a '.'
            The mitogenome sections assume that the files are named as above, $sample_id.$sequencing_technology.$date with
            additional information added to this with each proccess. The assembly process will add a .getorganelle${version}
            after the inital prefix and befor the file extension. Then the annotation process will add a .emma${version} to 
            the previous prefix.
        
        The 'else' statement is then just creating and empty channel if neither of the previous steps worked.
    */

    
    // Per-sample reference GenBanks (findMitoReference, from the assembly stage)
    // reused by the anthozoan annotation fixer. Partial / empty on precomputed or
    // skip paths; the annotation subworkflow re-resolves any sample missing one.
    ch_mitogenome_getorg_reference_gb = Channel.empty()
    ch_mitogenome_hifi_reference_gb   = Channel.empty()

    //
    // SUBWORKFLOW: MITOGENOME_ASSEMBLY_GETORG
    //

    if (!params.skip_mitogenome_assembly_getorg) {
        MITOGENOME_ASSEMBLY_GETORG (
            getorg_input,
            organelle_type // params.organelle_type
        )
        ch_mitogenome_getorg_assembly_fasta = MITOGENOME_ASSEMBLY_GETORG.out.assembly_fasta
        ch_mitogenome_getorg_assembly_log = MITOGENOME_ASSEMBLY_GETORG.out.assembly_log
        // Per-variant assembly results (first-pass + reseed + rgj) for the DB upload.
        ch_mitogenome_getorg_db_results = MITOGENOME_ASSEMBLY_GETORG.out.db_assembly_results
        ch_mitogenome_getorg_reference_gb = MITOGENOME_ASSEMBLY_GETORG.out.reference_gb
        ch_mitogenome_getorg_circularity_evidence = MITOGENOME_ASSEMBLY_GETORG.out.circularity_evidence
    } else if (params.precomputed_mitogenome_assembly_fasta_getorg) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_getorg_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_getorg, checkIfExists: false)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            def mt_assembly_prefix = parts.length > 2 ? parts[0..3].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        ch_mitogenome_getorg_assembly_log = Channel.fromPath(params.precomputed_mitogenome_assembly_log_getorg, checkIfExists: false)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            def mt_assembly_prefix = parts.length > 2 ? parts[0..3].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        // No check is re-run for precomputed assemblies, but the original run's
        // verdict is still on disk (*.getorg_check.tsv) -> reload it instead of
        // discarding it, so meta.circular doesn't silently regress to unknown.
        ch_mitogenome_getorg_circularity_evidence = Channel.fromPath(params.precomputed_mitogenome_circularity_evidence_getorg, checkIfExists: false)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            def mt_assembly_prefix = parts.length > 2 ? parts[0..3].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        // Precomputed inputs expose only the final assembly, so no per-variant
        // split is possible: fall back to a single DB row from the final fasta+log.
        ch_mitogenome_getorg_db_results = ch_mitogenome_getorg_assembly_fasta.join(ch_mitogenome_getorg_assembly_log, by: 0)
    } else {

        ch_mitogenome_getorg_assembly_fasta = Channel.empty()
        ch_mitogenome_getorg_assembly_log = Channel.empty()
        ch_mitogenome_getorg_circularity_evidence = Channel.empty()
        ch_mitogenome_getorg_db_results = Channel.empty()
    }

    if (!params.skip_mitogenome_assembly_getorg) {
        ch_assembly_summary_files = ch_assembly_summary_files.mix(MITOGENOME_ASSEMBLY_GETORG.out.summary_files)
    } else {
        ch_assembly_summary_files = ch_assembly_summary_files.mix(ch_mitogenome_getorg_assembly_fasta.map { meta, fasta -> fasta })
        ch_assembly_summary_files = ch_assembly_summary_files.mix(ch_mitogenome_getorg_assembly_log.map { meta, log -> log })
    }


    // 
    // SUBWORKFLOW: MITOGENOME_ASSEMBLY_MITOHIFI
    //

    if (!params.skip_mitogenome_assembly_hifi) {
        MITOGENOME_ASSEMBLY_MITOHIFI (
            mitohifi_input,
        )
        ch_mitogenome_hifi_assembly_fasta = MITOGENOME_ASSEMBLY_MITOHIFI.out.assembly_fasta
        ch_mitogenome_hifi_assembly_log = MITOGENOME_ASSEMBLY_MITOHIFI.out.assembly_log
        ch_mitogenome_hifi_reference_gb = MITOGENOME_ASSEMBLY_MITOHIFI.out.reference_gb
        ch_mitogenome_hifi_circularity_evidence = MITOGENOME_ASSEMBLY_MITOHIFI.out.circularity_evidence
        // Reference-free Oatk fallback contigs (empty channel unless the fallback is
        // enabled); folded into the annotation input below alongside the assemblies.
        ch_mitogenome_hifi_oatk_fasta = MITOGENOME_ASSEMBLY_MITOHIFI.out.oatk_fasta
        ch_mitogenome_hifi_oatk_log = MITOGENOME_ASSEMBLY_MITOHIFI.out.oatk_log
    } else if (params.precomputed_mitogenome_assembly_fasta_hifi) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_hifi_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_hifi, checkIfExists: false)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            def mt_assembly_prefix = parts.length > 2 ? parts[0..3].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        ch_mitogenome_hifi_assembly_log = Channel.fromPath(params.precomputed_mitogenome_assembly_log_hifi, checkIfExists: false)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            def mt_assembly_prefix = parts.length > 2 ? parts[0..3].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        // No circularity check is re-run for precomputed assemblies, but the
        // original run's verdict is still on disk (*.circularity_check.tsv) ->
        // reload it instead of discarding it, so meta.circular doesn't silently
        // regress to unknown.
        ch_mitogenome_hifi_circularity_evidence = Channel.fromPath(params.precomputed_mitogenome_circularity_evidence_hifi, checkIfExists: false)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            def mt_assembly_prefix = parts.length > 2 ? parts[0..3].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        ch_mitogenome_hifi_oatk_fasta = Channel.empty()
        ch_mitogenome_hifi_oatk_log = Channel.empty()
    } else {

        ch_mitogenome_hifi_assembly_fasta = Channel.empty()
        ch_mitogenome_hifi_assembly_log = Channel.empty()
        ch_mitogenome_hifi_circularity_evidence = Channel.empty()
        ch_mitogenome_hifi_oatk_fasta = Channel.empty()
        ch_mitogenome_hifi_oatk_log = Channel.empty()
    }

    if (!params.skip_mitogenome_assembly_hifi) {
        ch_assembly_summary_files = ch_assembly_summary_files.mix(MITOGENOME_ASSEMBLY_MITOHIFI.out.summary_files)
    } else {
        ch_assembly_summary_files = ch_assembly_summary_files.mix(ch_mitogenome_hifi_assembly_fasta.map { meta, fasta -> fasta })
        ch_assembly_summary_files = ch_assembly_summary_files.mix(ch_mitogenome_hifi_assembly_log.map { meta, log -> log })
    }

    //
    // mix the hifi and getorganelle outputs to feed into the annotation
    //
    // Samples whose assembler finished without producing a contig carry an
    // empty placeholder FASTA. Filter those out so annotation / LCA / QC only
    // run on real assemblies; PUSH_MTDNA_ASSM_RESULTS below still records the
    // failure in the database.
    ch_all_assembly_fasta = ch_mitogenome_hifi_assembly_fasta
        .mix(ch_mitogenome_getorg_assembly_fasta)
        .mix(ch_mitogenome_hifi_oatk_fasta)
    ch_failed_assembly_fasta = ch_all_assembly_fasta
        .filter { _meta, fasta -> fasta.size() == 0 }
    ch_annotation_input = ch_all_assembly_fasta
        .filter { _meta, fasta -> fasta.size() > 0 }

    // Auto-curation: collapse a clean head-to-tail concatemer (an assembly ~2x
    // the true length, detected by the circularity check) to a single monomer
    // before annotation. collapse_concatemer.py re-verifies the multimer by
    // self-alignment and passes anything it cannot confirm through unchanged, so
    // this can only ever fix a genuine duplication, never corrupt an assembly.
    // Assemblies with no circularity evidence (remainder) bypass the step.
    ch_collapse_branched = ch_annotation_input
        .join(
            ch_mitogenome_hifi_circularity_evidence.mix(ch_mitogenome_getorg_circularity_evidence),
            by: 0, remainder: true
        )
        .branch { _meta, fasta, evidence ->
            concatemer: evidence != null && fasta != null && evidenceAnomalyType(evidence) == 'concatemer'
            bypass: true
        }

    COLLAPSE_CONCATEMER(
        ch_collapse_branched.concatemer.map { meta, fasta, evidence -> [meta, fasta, evidence] }
    )

    // The collapse report feeds the assembly summary so a collapsed concatemer is
    // reported at its monomer length rather than re-flagged as over-length.
    ch_assembly_summary_files = ch_assembly_summary_files.mix(
        COLLAPSE_CONCATEMER.out.report.map { _meta, report -> report }
    )
    ch_versions = ch_versions.mix(COLLAPSE_CONCATEMER.out.versions.first())

    // Rebuild the annotation input from the (possibly collapsed) FASTAs plus the
    // assemblies that had no evidence to collapse against.
    ch_collapsed_canonical_fasta = COLLAPSE_CONCATEMER.out.fasta
        .join(COLLAPSE_CONCATEMER.out.evidence, by: 0)
        .map { meta, fasta, evidence ->
            def circular = evidenceFinalVerdictCircular(evidence)
            [ circular == null ? meta : meta + [ circular: circular ], fasta ]
        }
    ch_canonical_assembly_fasta = ch_collapsed_canonical_fasta
        .mix(ch_collapse_branched.bypass.map { meta, fasta, _evidence -> [meta, fasta] })
        .filter { _meta, fasta -> fasta != null && fasta.size() > 0 }

    ch_canonical_circularity_evidence = COLLAPSE_CONCATEMER.out.evidence
        .mix(ch_collapse_branched.bypass
            .filter { _meta, _fasta, evidence -> evidence != null }
            .map { meta, _fasta, evidence -> [meta, evidence] })

    // Assemblies below the configured biological minimum are retained for SQL
    // and summary reporting, but do not enter EMMA/MITOS/table2asn.
    ch_annotation_input = ch_canonical_assembly_fasta
        .filter { _meta, fasta -> fastaSequenceLength(fasta) >= params.mitogenome_summary_min_length }

    // Sanitise FASTA before annotation to avoid duplicate IDs / multi-contig issues
    SANITISE_FASTA(
        ch_annotation_input
    )
    ch_annotation_input_sanitised = SANITISE_FASTA.out

    //
    // SUBWORKFLOW: MITOGENOME_ANNOTATION
    //
    // Per-sample reference GenBanks from whichever assembler ran, for the
    // anthozoan annotation fixer to reuse before re-downloading.
    ch_annotation_reference_gb = ch_mitogenome_hifi_reference_gb
        .mix(ch_mitogenome_getorg_reference_gb)

    if (!params.skip_mitogenome_annotation) {
        MITOGENOME_ANNOTATION (
            ch_annotation_input_sanitised,
            curated_blast_db,
            nt_blast_db,
            mitos_refdb,
            ch_annotation_reference_gb
        )
        ch_mitogenome_annotation_results = MITOGENOME_ANNOTATION.out.annotation_results
        ch_mitogenome_blast_results = MITOGENOME_ANNOTATION.out.blast_filtered_results
        ch_mitogenome_lca_results = MITOGENOME_ANNOTATION.out.lca_results
        ch_mitogenome_lca_raw_results = MITOGENOME_ANNOTATION.out.lca_raw_results

        // Feed the per-sample reference-relevance flag into the assembly summary so
        // a wrong-family reference surfaces as a manual_review_reason.
        ch_assembly_summary_files = ch_assembly_summary_files.mix(
            MITOGENOME_ANNOTATION.out.reference_relevance.map { meta, f -> f })
    } else if (params.precomputed_mitogenome_annotation_results) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_annotation_results = Channel.fromPath(params.precomputed_mitogenome_annotation_results)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def sequencing_type = parts.length > 1 ? parts[1] : null
            def date = parts.length > 2 ? parts[2] : null
            // Reconstruct the assembly prefix (id.tech.date.assembler) from the
            // annotation filename so downstream naming (e.g. <prefix>.annotation_stats.csv)
            // stays per-assembly. Without this, meta.mt_assembly_prefix is null
            // on the upload-only path and CSVs collide as null.annotation_stats.csv.
            def mt_assembly_prefix = parts.length > 3 ? parts[0..3].join('.') : filename
            return [ mt_assembly_prefix, [meta_id, sequencing_type, date], file ]
        }
        // FORMAT_FILES needs the whole per-assembly annotation bundle (fasta, gff,
        // tbl/gb) staged together -- glob matches every one of those files
        // individually, so group them back into one list per assembly prefix,
        // mirroring the path("annotation/*") bundle that
        // MITOGENOME_ANNOTATION.out.annotation_results emits when annotation isn't
        // skipped. Without this, only the last-matched file (e.g. the GFF) reaches
        // FORMAT_FILES and it fails with "Missing in .: FASTA".
        .groupTuple(by: 0)
        .map { mt_assembly_prefix, sample_keys, files -> [ sample_keys[0], mt_assembly_prefix, files ] }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, files, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, files)
        }
        ch_mitogenome_blast_results = Channel.fromPath(params.precomputed_mitogenome_blast_results)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            // blast.<gene_type>.<og_id>.<tech>.<date>.<assembler>[.<annot_version>].filtered
            def meta_id = parts[2]
            def sequencing_type = parts.length > 1 ? parts[3] : null
            def date = parts.length > 2 ? parts[4] : null
            // Reconstruct the assembly prefix (id.tech.date.assembler) so downstream
            // naming (e.g. <prefix>.lca_blast.upload.txt) stays per-assembly. Without
            // this, meta.mt_assembly_prefix is null on the upload-only path and
            // multiple assembly attempts for one OG collide on the same filename.
            def mt_assembly_prefix = parts.length > 5 ? parts[2..5].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        ch_mitogenome_lca_results = Channel.fromPath(params.precomputed_mitogenome_lca_results)
        .map { file ->
            def filename = file.baseName
            def parts = filename.split('\\.')
            // lca.<gene_type>.<og_id>.<tech>.<date>.<assembler>[.<annot_version>]
            def meta_id = parts[2]
            def sequencing_type = parts.length > 3 ? parts[3] : null
            def date = parts.length > 4 ? parts[4] : null
            def mt_assembly_prefix = parts.length > 5 ? parts[2..5].join('.') : filename
            return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
        }
        .combine(ch_samplesheet_meta, by: 0)
        .map { sample_key, mt_assembly_prefix, file, meta ->
            def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
            return tuple(meta_ext, file)
        }
        ch_mitogenome_lca_raw_results = params.precomputed_mitogenome_lca_raw_results
            ? Channel.fromPath(params.precomputed_mitogenome_lca_raw_results)
                .map { file ->
                    def filename = file.baseName
                    def parts = filename.split('\\.')
                    // lca_raw.<region>.<og_id>.<tech>.<date>.<assembler>[.<annot_version>]
                    def meta_id = parts[2]
                    def sequencing_type = parts.length > 3 ? parts[3] : null
                    def date = parts.length > 4 ? parts[4] : null
                    def mt_assembly_prefix = parts.length > 5 ? parts[2..5].join('.') : filename
                    return [ [meta_id, sequencing_type, date], mt_assembly_prefix, file ]
                }
                .combine(ch_samplesheet_meta, by: 0)
                .map { sample_key, mt_assembly_prefix, file, meta ->
                    def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
                    return tuple(meta_ext, file)
                }
            : Channel.empty()
    } else {

        ch_mitogenome_annotation_results = Channel.empty()
        ch_mitogenome_blast_results = Channel.empty()
        ch_mitogenome_lca_results = Channel.empty()
        ch_mitogenome_lca_raw_results = Channel.empty()
    }

    //
    // Combine outputs for data uploads
    //add in mitohifi stuff too

    // Preserve every GetOrganelle variant in SQL, while replacing the selected
    // final variant (same mt_assembly_prefix) with its canonical post-curation
    // FASTA/meta. MitoHiFi/Oatk canonical results have no raw GetOrg variant and
    // are therefore added normally. Grouping by prefix prevents duplicate rows.
    ch_canonical_assembly_logs = ch_mitogenome_getorg_assembly_log
        .mix(ch_mitogenome_hifi_assembly_log)
        .mix(ch_mitogenome_hifi_oatk_log)
        .map { meta, log -> [ meta.mt_assembly_prefix, log ] }

    ch_canonical_upload_candidates = ch_canonical_assembly_fasta
        .mix(ch_failed_assembly_fasta)
        .map { meta, fasta -> [ meta.mt_assembly_prefix, meta, fasta ] }
        .join(ch_canonical_assembly_logs, by: 0)
        .map { prefix, meta, fasta, log ->
            [ prefix, [ priority: 1, meta: meta, fasta: fasta, log: log ] ]
        }

    ch_raw_getorg_upload_candidates = ch_mitogenome_getorg_db_results
        .map { meta, fasta, log ->
            [ meta.mt_assembly_prefix, [ priority: 0, meta: meta, fasta: fasta, log: log ] ]
        }

    ch_mitogenome_assembly_results = ch_raw_getorg_upload_candidates
        .mix(ch_canonical_upload_candidates)
        .groupTuple(by: 0)
        .map { _prefix, candidates ->
            def selected = candidates.max { it.priority }
            [ selected.meta, selected.fasta, selected.log ]
        }

    //
    // SUBWORKFLOW: UPLOAD_RESULTS
    //
    // Conditional uploading of results to SQL and species check - only run if not skipped
    // All these processes access the OceanOmics PostgreSQL database.
    def ch_qc_input = Channel.empty()
    if (!params.skip_upload_results && params.sql_config) {
        // Per-sample circularity-check evidence from both assemblers feeds the QC
        // gate (anomaly block; the circular condition itself comes via meta.circular).
        ch_mitogenome_circularity_evidence = ch_canonical_circularity_evidence

        UPLOAD_RESULTS (
            ch_mitogenome_assembly_results,
            ch_mitogenome_annotation_results,
            ch_mitogenome_blast_results,
            ch_mitogenome_lca_results,
            ch_mitogenome_lca_raw_results,
            ch_mitogenome_circularity_evidence,
            sql_config // params.sql_config
        )

        ch_qc_input = UPLOAD_RESULTS.out.qc_ready.join(ch_mitogenome_annotation_results, by:0)

        // If the LCA validation is correct, then run the QC to prepare for submission to GenBank
        // Need to add this into the pipeline.
        // Now that protein lengths are being added to the database it could provide a list of 
        // non submitted mitogenomes they can be grouped with to submit and then say when there is a 
        // group of similar mitogenomes they can be submitted as a batch.
        MITOGENOME_QC (
            ch_qc_input // tuple val(meta), val(species_name), val(proceed_qc true/false), val(circular true/false), path(annotation/*)
        )
        UPLOAD_ENA_RESULTS (
            MITOGENOME_QC.out.ena_validation_records,
            UPLOAD_RESULTS.out.upload_status_files,
            sql_config
        )
        ch_assembly_summary_files = ch_assembly_summary_files.mix(UPLOAD_RESULTS.out.assembly_summary_files)
    } else if (!params.skip_upload_results && !params.sql_config) {
        log.warn "Skipping upload/QC because --sql_config not provided"
    }

    //
    // MODULE: Mitogenome assembly summary for MultiQC custom content
    //

    ch_assembly_summary_inputs = ch_assembly_summary_files.ifEmpty(
        file("$projectDir/assets/multiqc_config.yml", checkIfExists: false)
    )

    MITOGENOME_ASSEMBLY_SUMMARY (
        ch_assembly_summary_inputs.collect()
    )

    ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ASSEMBLY_SUMMARY.out.table)
    ch_versions = ch_versions.mix(MITOGENOME_ASSEMBLY_SUMMARY.out.versions)

    //
    // Collect all MultiQC files from all subworkflows
    //

    if (!params.skip_mitogenome_assembly_getorg) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ASSEMBLY_GETORG.out.multiqc_files)}
    if (!params.skip_mitogenome_assembly_hifi) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ASSEMBLY_MITOHIFI.out.multiqc_files)}
    if (!params.skip_mitogenome_annotation) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ANNOTATION.out.multiqc_files)}
    if (!params.skip_upload_results && params.sql_config) {ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_RESULTS.out.multiqc_files)}
    if (!params.skip_upload_results && params.sql_config) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_QC.out.multiqc_files)}
    if (!params.skip_upload_results && params.sql_config) {ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_ENA_RESULTS.out.multiqc_files)}

    // 
    // Collect all versions from subworkflows
    //

    if (!params.skip_mitogenome_assembly_getorg) {
        ch_versions = ch_versions.mix(
            MITOGENOME_ASSEMBLY_GETORG.out.versions,
        )
    }
    if (!params.skip_mitogenome_assembly_hifi) {
        ch_versions = ch_versions.mix(
            MITOGENOME_ASSEMBLY_MITOHIFI.out.versions
        )
    }
    if (!params.skip_mitogenome_annotation) {
        ch_versions = ch_versions.mix(MITOGENOME_ANNOTATION.out.versions)
    }
    // Upload + QC subworkflows provide versions; guard independently
    if (!params.skip_upload_results && params.sql_config) {
        ch_versions = ch_versions.mix(UPLOAD_RESULTS.out.versions)
    }
    // Run of MITOGENOME_QC depends on upstream evaluation; include if present
    if (!params.skip_upload_results && params.sql_config) {
        ch_versions = ch_versions.mix(MITOGENOME_QC.out.versions)
    }
    if (!params.skip_upload_results && params.sql_config) {
        ch_versions = ch_versions.mix(UPLOAD_ENA_RESULTS.out.versions)
    }







    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oceangenomesmitogenomes_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: false)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: false) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: false) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: false) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: false)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files_collected = ch_multiqc_files.collect()
    ch_multiqc_config_list = ch_multiqc_config.toList()
    ch_multiqc_custom_config_list = ch_multiqc_custom_config.toList()
    ch_multiqc_logo_list = ch_multiqc_logo.toList()

    if (!params.skip_per_sample_multiqc) {
        MULTIQC_PER_SAMPLE (
            ch_multiqc_files_collected,
            ch_multiqc_config_list,
            ch_multiqc_custom_config_list,
            ch_multiqc_logo_list,
            [],
            []
        )
    }

    MULTIQC (
        ch_multiqc_files_collected,
        ch_multiqc_config_list,
        ch_multiqc_custom_config_list,
        ch_multiqc_logo_list,
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
