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
include { MIRROR_MTDNA_TO_COLLAPSED } from '../modules/local/mirror_mtdna_to_collapsed'
include { UPLOAD_RESULTS; UPLOAD_ENA_RESULTS } from '../subworkflows/local/upload_results_mito'
include { MITOGENOME_QC             } from '../subworkflows/local/mitogenome_qc'
include { SANITISE_FASTA           } from '../modules/local/sanitise_fasta/main'
include { MITOGENOME_COVERAGE      } from '../modules/local/mitogenome_coverage/main'
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

// Build the per-assembly SQL upload rows -- every canonical (post-curation) assembly plus
// every GetOrganelle provenance variant -- each carrying its uniform depth TSV, or the
// header-only placeholder where none was measured. Emits each row as soon as THAT assembly
// is done, never waiting for the rest of the run.
//
// canonical_rows / variant_rows are [ mt_assembly_prefix, meta, fasta, log ];
// mito_depth is [ mt_assembly_prefix, depth.tsv ]. Returns [ meta, fasta, log, depth ].
//
// Extracted from the workflow body so tests/assembly_upload_streaming can exercise the real
// implementation rather than a copy of it -- see the call site for the run this protects.
def buildAssemblyUploadRows(canonical_rows, variant_rows, mito_depth, no_depth_file, min_length) {
    // Rows that are never measured BY CONSTRUCTION: an empty (failed) assembly and an
    // under-length one never reach annotation, and a provenance variant is not the molecule
    // that was annotated. Both predicates are row-local -- the same ones used to build
    // ch_annotation_input -- so routing them around the depth join costs nothing. It matters
    // because a remainder join can only classify a row as unmatched once its inputs CLOSE,
    // which would pin every failed and every provenance row's SQL upload to end of run.
    def routed = canonical_rows.branch { _prefix, _meta, fasta, _log ->
        measured:   fasta.size() > 0 && fastaSequenceLength(fasta) >= min_length
        unmeasured: true
    }

    def unmeasured = routed.unmeasured
        .mix(variant_rows)
        .map { _prefix, meta, fasta, log -> [ meta, fasta, log, no_depth_file ] }

    // The remaining arm keeps remainder:true for one case only: a mixed run where an
    // assembler is skipped but its results are precomputed has no reads channel for those
    // samples, so MITOGENOME_COVERAGE never runs on them and they legitimately have no
    // depth. On the live path every measurable assembly has reads, so the join matches and
    // emits immediately -- and those are the rows the QC gate waits on.
    def measured = routed.measured
        .join(mito_depth, by: 0, remainder: true)
        .filter { items -> items[1] != null }   // keep assembly rows; drop depth-only remainder
        .map { items ->
            def depth = (items.size() > 4 && items[4] != null) ? items[4] : no_depth_file
            [ items[1], items[2], items[3], depth ]
        }

    return measured.mix(unmeasured)
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

// Read the 'action' column from a COLLAPSE_CONCATEMER report (concatemer_collapse.tsv).
// 'collapsed' means the assembly was genuinely rewritten to its monomer; anything
// else (e.g. 'passthrough') means it was left unchanged.
def collapseAction(tsv) {
    try {
        def rows = tsv.text.readLines()
        if (rows.size() < 2) return ''
        def header = rows[0].split('\t', -1)
        def idx = header.findIndexOf { it.trim() == 'action' }
        def values = rows[1].split('\t', -1)
        return idx >= 0 && idx < values.size() ? values[idx].trim().toLowerCase() : ''
    } catch (ignored) {
        return ''
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

    // Stand-in for an assembly with no circularity check on disk. The assembly subworkflows
    // guarantee one evidence row per assembly they emit; the precomputed paths below have to
    // reconstruct that guarantee from whatever *_check.tsv files happen to be present.
    def no_circularity_evidence = file("${projectDir}/assets/empty_circularity_check.tsv", checkIfExists: true)

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

    // Per-sample bundles of the assembly-stage mtdna files (keyed by assembly
    // prefix), used to mirror the whole mtdna folder into <prefix>_collapsed/mtdna
    // for genuinely collapsed samples. Empty on precomputed / skip paths (nothing
    // was assembled this run to mirror).
    ch_mitogenome_getorg_mtdna_files = Channel.empty()
    ch_mitogenome_hifi_mtdna_files   = Channel.empty()

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
        ch_mitogenome_getorg_mtdna_files = MITOGENOME_ASSEMBLY_GETORG.out.mtdna_files
    } else if (params.precomputed_mitogenome_assembly_fasta_getorg) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_getorg_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_getorg, checkIfExists: false)
        // Keep each assembly FASTA only from the variant directory it belongs to
        // (<prefix>/mtdna/<prefix>.fasta). A genuinely collapsed sample also publishes
        // its <prefix>_collapsed.fasta back into the original <prefix>/mtdna dir for
        // provenance, so the greedy re-glob would otherwise pick the same basename up
        // from two mtdna dirs and collide in MITOGENOME_ASSEMBLY_SUMMARY's flat stageAs.
        .filter { it.baseName == it.parent.parent.name }
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
        // MIRROR_MTDNA_TO_COLLAPSED restages a collapsed sample's check into
        // <prefix>_collapsed/mtdna, so the glob matches the same evidence twice for those
        // samples. Keep only the copy under its own <prefix>/mtdna dir: the evidence is now
        // joined one-to-one against the assembly, so a duplicate key would annotate the
        // sample twice. Same directory-shape filter the oatk fasta glob below uses.
        .filter { f -> f.parent.parent.name == f.baseName.replaceAll(/\.[^.]+$/, '') }
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
        // Restore the one-row-per-assembly contract the assembly subworkflows guarantee:
        // a precomputed run may have no *_check.tsv for some assemblies (the glob uses
        // checkIfExists: false), and downstream now attaches evidence with a plain per-key
        // join, so an absent row would drop the sample at the QC gate rather than delay it.
        // remainder: true is safe HERE, unlike the collapse join below, because both sides are
        // Channel.fromPath globs that close immediately -- there is no running task to wait on.
        ch_mitogenome_getorg_circularity_evidence = ch_mitogenome_getorg_assembly_fasta
            .map { meta, _fasta -> [ meta.mt_assembly_prefix, meta ] }
            .join(
                ch_mitogenome_getorg_circularity_evidence.map { meta, ev -> [ meta.mt_assembly_prefix, ev ] },
                by: 0, remainder: true
            )
            .filter { items -> items[1] != null }   // keep assembly rows; drop evidence-only remainder
            .map { items -> [ items[1], (items.size() > 2 && items[2] != null) ? items[2] : no_circularity_evidence ] }

        // Precomputed inputs expose only the final assembly, so there are no provenance
        // variants to record: the one row a precomputed sample could produce carries the
        // same mt_assembly_prefix as its canonical row below, which supersedes it. This
        // channel therefore contributes nothing on this path -- it was already discarded
        // by the selection that used to run in the parent, and the plain mix that replaced
        // that selection has no way to discard it, so leave it empty.
        ch_mitogenome_getorg_db_results = Channel.empty()
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
        ch_mitogenome_hifi_mtdna_files = MITOGENOME_ASSEMBLY_MITOHIFI.out.mtdna_files
        // Reference-free Oatk fallback contigs (empty channel unless the fallback is
        // enabled); folded into the annotation input below alongside the assemblies.
        ch_mitogenome_hifi_oatk_fasta = MITOGENOME_ASSEMBLY_MITOHIFI.out.oatk_fasta
        ch_mitogenome_hifi_oatk_log = MITOGENOME_ASSEMBLY_MITOHIFI.out.oatk_log
    } else if (params.precomputed_mitogenome_assembly_fasta_hifi) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_hifi_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_hifi, checkIfExists: false)
        // Keep each assembly FASTA only from the variant directory it belongs to
        // (<prefix>/mtdna/<prefix>.fasta). A genuinely collapsed sample also publishes
        // its <prefix>_collapsed.fasta back into the original <prefix>/mtdna dir for
        // provenance, so the greedy re-glob would otherwise pick the same basename up
        // from two mtdna dirs and collide in MITOGENOME_ASSEMBLY_SUMMARY's flat stageAs.
        // Also covers the oatk fallback FASTAs, which ride this same hifi precomputed path.
        .filter { it.baseName == it.parent.parent.name }
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
        // Drop the <prefix>_collapsed/mtdna mirror copy -- see the GetOrganelle glob above.
        .filter { f -> f.parent.parent.name == f.baseName.replaceAll(/\.[^.]+$/, '') }
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
        // Oatk fallback assemblies are published beside the MitoHiFi results, so a
        // skipped/precomputed run must reload them too, otherwise every oatk-rescued
        // sample silently drops out of annotation, the canonical logs and the upload.
        // Empty (no matches) when the fallback never ran. Same variant-dir filter as
        // the hifi fasta above: the *oatk*.fasta glob also matches the native
        // <prefix>.mito.ctg.fasta provenance copy, which this excludes so only the
        // canonical <prefix>.fasta rides through.
        ch_mitogenome_hifi_oatk_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_oatk, checkIfExists: false)
        .filter { it.baseName == it.parent.parent.name }
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
        ch_mitogenome_hifi_oatk_log = Channel.fromPath(params.precomputed_mitogenome_assembly_log_oatk, checkIfExists: false)
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

        // Restore the one-row-per-assembly contract the assembly subworkflows guarantee, for
        // the MitoHiFi and oatk assemblies alike -- both feed ch_all_assembly_fasta below, so
        // both need evidence. See the GetOrganelle equivalent above for why remainder: true is
        // safe on the precomputed path and not downstream.
        ch_mitogenome_hifi_circularity_evidence = ch_mitogenome_hifi_assembly_fasta
            .mix(ch_mitogenome_hifi_oatk_fasta)
            .map { meta, _fasta -> [ meta.mt_assembly_prefix, meta ] }
            .join(
                ch_mitogenome_hifi_circularity_evidence.map { meta, ev -> [ meta.mt_assembly_prefix, ev ] },
                by: 0, remainder: true
            )
            .filter { items -> items[1] != null }   // keep assembly rows; drop evidence-only remainder
            .map { items -> [ items[1], (items.size() > 2 && items[2] != null) ? items[2] : no_circularity_evidence ] }
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
        // Oatk fallbacks are part of the live subworkflow's summary_files, so include
        // the reloaded oatk fasta + log here too, or skipped runs would summarise every
        // sample except the oatk-rescued ones.
        ch_assembly_summary_files = ch_assembly_summary_files.mix(ch_mitogenome_hifi_oatk_fasta.map { meta, fasta -> fasta })
        ch_assembly_summary_files = ch_assembly_summary_files.mix(ch_mitogenome_hifi_oatk_log.map { meta, log -> log })
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
    // Assemblies whose evidence records no concatemer bypass the step.
    //
    // A plain join, not join(..., remainder: true). Both assembly subworkflows now guarantee
    // one evidence row per assembly they emit (empty_circularity_check.tsv where no check
    // ran), so there is nothing left for the remainder path to catch. That matters because a
    // remainder join can only classify an UNMATCHED row once the evidence channel closes:
    // every sample without evidence was held back until the slowest assembly in the run
    // finished, reaching annotation as an end-of-run batch instead of incrementally.
    //
    // Keyed on mt_assembly_prefix rather than the whole meta. Assembly FASTA and evidence
    // reach here from different processes and, on a -resume, from different cache entries,
    // so a restored meta may no longer `equals` the live one -- the failure already
    // documented for the oatk reference join, and the reason OG109/OG2089 once never reached
    // OATK_CHECK. Under the old remainder join a key miss merely delayed a sample down the
    // bypass path; under a plain join it would drop it from annotation entirely, so the
    // robust key is now load-bearing. The prefix is unique per assembly (the duplicate
    // <prefix>_collapsed publish copy exists only on disk, never in this channel).
    ch_collapse_branched = ch_annotation_input
        .map { meta, fasta -> [ meta.mt_assembly_prefix, meta, fasta ] }
        .join(
            ch_mitogenome_hifi_circularity_evidence
                .mix(ch_mitogenome_getorg_circularity_evidence)
                .map { meta, evidence -> [ meta.mt_assembly_prefix, evidence ] },
            by: 0
        )
        .map { _prefix, meta, fasta, evidence -> [ meta, fasta, evidence ] }
        .branch { _meta, _fasta, evidence ->
            concatemer: evidenceAnomalyType(evidence) == 'concatemer'
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

    // Provenance mirror. For a GENUINELY collapsed sample (action == collapsed in
    // the report), restage the whole assembly-stage mtdna folder (original assembly
    // + circularity check + collapse artefacts) into the <prefix>_collapsed/mtdna
    // dir that downstream annotation forks into (via the _collapsed FASTA basename),
    // so the curated variant directory carries the full provenance in one place.
    // Passthrough samples are untouched and keep only their original mtdna dir.
    // Keyed on the original assembly prefix (mt_assembly_prefix), which is stable
    // across the collapse -- the _collapsed suffix lives only in the FASTA basename.
    ch_collapse_bundle = COLLAPSE_CONCATEMER.out.fasta
        .join(COLLAPSE_CONCATEMER.out.report, by: 0)
        .join(COLLAPSE_CONCATEMER.out.evidence, by: 0)
        .filter { _meta, _fasta, report, _evidence -> collapseAction(report) == 'collapsed' }
        .map { meta, fasta, report, evidence -> [ meta.mt_assembly_prefix, meta, fasta, report, evidence ] }

    ch_mtdna_bundle = ch_mitogenome_hifi_mtdna_files.mix(ch_mitogenome_getorg_mtdna_files)

    ch_mirror_input = ch_collapse_bundle
        .join(ch_mtdna_bundle, by: 0, remainder: true)
        .filter { items -> items[1] != null }   // keep collapse rows; drop bundle-only remainder
        .map { items ->
            def meta      = items[1]
            def fasta     = items[2]   // <prefix>_collapsed.fasta (renamed by COLLAPSE_CONCATEMER)
            def report    = items[3]
            def evidence  = items[4]
            def bundle    = (items.size() > 5 && items[5] != null) ? items[5] : []
            def newPrefix = fasta.baseName   // <prefix>_collapsed -> drives the mirror's publish dir
            [ meta + [ mt_assembly_prefix: newPrefix ], [ fasta, report, evidence ] + bundle ]
        }

    MIRROR_MTDNA_TO_COLLAPSED( ch_mirror_input )

    // Rebuild the annotation input from the (possibly collapsed) FASTAs plus the
    // assemblies that had no concatemer to collapse.
    ch_collapsed_canonical_fasta = COLLAPSE_CONCATEMER.out.fasta
        .join(COLLAPSE_CONCATEMER.out.evidence, by: 0)
        .map { meta, fasta, evidence ->
            def circular = evidenceFinalVerdictCircular(evidence)
            [ circular == null ? meta : meta + [ circular: circular ], fasta ]
        }
    ch_canonical_assembly_fasta = ch_collapsed_canonical_fasta
        .mix(ch_collapse_branched.bypass.map { meta, fasta, _evidence -> [meta, fasta] })
        .filter { _meta, fasta -> fasta.size() > 0 }

    // One evidence row per canonical assembly, keyed by mt_assembly_prefix for the QC gate.
    //
    // Prefix-keyed, not meta-keyed: the collapsed branch above rewrites meta (it folds the
    // post-collapse verdict into meta.circular) while COLLAPSE_CONCATEMER.out.evidence still
    // carries the pre-collapse meta, so a whole-meta lookup misses every genuinely collapsed
    // sample and silently hands it the "no anomaly" stand-in. mt_assembly_prefix is stable
    // across the collapse -- the _collapsed suffix lives only in the FASTA basename.
    ch_canonical_circularity_evidence = COLLAPSE_CONCATEMER.out.evidence
        .mix(ch_collapse_branched.bypass.map { meta, _fasta, evidence -> [meta, evidence] })
        .map { meta, evidence -> [ meta.mt_assembly_prefix, evidence ] }

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
    // MODULE: MITOGENOME_COVERAGE -- one uniform, cross-platform depth number.
    //
    //   Remaps each sample's own full read set to the assembly and reports mean
    //   per-base depth. This replaces three quantities that were never comparable:
    //   GetOrganelle's k-mer coverage (~0.2x true depth, off a reduced read set),
    //   MitoHiFi's depth over only the reads a related-species reference recruited
    //   (so a divergent reference depressed it), and Oatk's absence of any number.
    //
    //   Deliberately placed on the SANITISE_FASTA output, i.e. the exact molecule
    //   that reaches annotation and GenBank: post-collapse (measuring a concatemer
    //   would halve the depth), post-reseed (only the variant that won matters), and
    //   post-concatenation (SANITISE_FASTA rewrites a multi-contig assembly into a
    //   single _concat record). Anything that failed or fell below the length floor
    //   was already filtered out above, so no remap is spent on a dead assembly.
    //
    ch_depth_reads = Channel.empty()
    if (!params.skip_mitogenome_assembly_getorg) {
        ch_depth_reads = ch_depth_reads.mix(MITOGENOME_ASSEMBLY_GETORG.out.depth_reads)
    }
    if (!params.skip_mitogenome_assembly_hifi) {
        ch_depth_reads = ch_depth_reads.mix(MITOGENOME_ASSEMBLY_MITOHIFI.out.depth_reads)
    }

    ch_mito_depth = Channel.empty()
    if (!params.skip_mitogenome_depth) {
        // Inner join on mt_assembly_prefix: a precomputed / assembly-skipped run has
        // no reads channel at all, and those samples simply get no depth (the SQL
        // push falls back to the empty placeholder rather than being dropped).
        MITOGENOME_COVERAGE (
            ch_annotation_input_sanitised
                .map { meta, fasta -> [ meta.mt_assembly_prefix, meta, fasta ] }
                .join(ch_depth_reads, by: 0)
                .map { _prefix, meta, fasta, reads -> [ meta, fasta, reads ] }
        )
        ch_mito_depth = MITOGENOME_COVERAGE.out.depth
        ch_assembly_summary_files = ch_assembly_summary_files.mix(
            ch_mito_depth.map { _meta, tsv -> tsv })
        ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_COVERAGE.out.tool_params.collect { it[1] })
        ch_versions = ch_versions.mix(MITOGENOME_COVERAGE.out.versions.first())
    }

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
        ch_mitogenome_region_counts = MITOGENOME_ANNOTATION.out.region_counts

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

        // Region count per sample for the downstream group sizing. On this path
        // the results are already on disk, so counting the precomputed BLAST
        // files per sample is exact and costs nothing: the fromPath channel
        // closes at startup, long before anything it gates.
        ch_mitogenome_region_counts = ch_mitogenome_blast_results
            .map { meta, _file -> [ meta, 1 ] }
            .groupTuple(by: 0)
            .map { meta, ones -> [ meta, ones.size() ] }
    } else {

        ch_mitogenome_annotation_results = Channel.empty()
        ch_mitogenome_blast_results = Channel.empty()
        ch_mitogenome_lca_results = Channel.empty()
        ch_mitogenome_lca_raw_results = Channel.empty()
        ch_mitogenome_region_counts = Channel.empty()
    }

    //
    // Combine outputs for data uploads
    //add in mitohifi stuff too

    // Preserve every GetOrganelle provenance variant in SQL alongside the canonical
    // post-curation row for each assembly. MitoHiFi/Oatk results have no GetOrganelle
    // variant and are added the same way.
    //
    // This used to mix the two channels and reconcile them with groupTuple(by: prefix),
    // picking the canonical row where a GetOrganelle first-pass variant collided with it.
    // A bare groupTuple emits NOTHING until its source channel closes, and that source runs
    // back to both assembly subworkflows -- so a single queued assembly task anywhere in the
    // run left PUSH_MTDNA_ASSM_RESULTS with zero tasks, which in turn emptied qc_ready (it
    // gates on the upload receipt) and stopped every finished sample from reaching QC. Run
    // mitogenomes-missing-audit-5 hit exactly that: 141 samples finished, one REFERENCE_RANK
    // sat behind a maintenance reservation, and nothing progressed past QC_SUMMARY.
    //
    // The fix is upstream disambiguation rather than a cleverer operator here: the
    // GetOrganelle subworkflow no longer emits the first-pass variant that collided (its
    // basename IS the sample's mt_assembly_prefix, so the canonical row always superseded it
    // anyway), and the precomputed path emits nothing at all. With no collisions left there
    // is nothing to reconcile, so a plain mix suffices and each row reaches SQL as soon as
    // its OWN assembly is done.
    ch_canonical_assembly_logs = ch_mitogenome_getorg_assembly_log
        .mix(ch_mitogenome_hifi_assembly_log)
        .mix(ch_mitogenome_hifi_oatk_log)
        .map { meta, log -> [ meta.mt_assembly_prefix, log ] }

    ch_canonical_upload_rows = ch_canonical_assembly_fasta
        .mix(ch_failed_assembly_fasta)
        .map { meta, fasta -> [ meta.mt_assembly_prefix, meta, fasta ] }
        .join(ch_canonical_assembly_logs, by: 0)

    ch_variant_upload_rows = ch_mitogenome_getorg_db_results
        .map { meta, fasta, log -> [ meta.mt_assembly_prefix, meta, fasta, log ] }

    // Attach the uniform depth TSV, keyed on the assembly prefix. Rows with no depth are
    // still uploaded, just with the header-only placeholder (push_mtdna_assm_results.py
    // records them as 'not_measured'). See buildAssemblyUploadRows for the routing.
    def no_depth_file = file("${projectDir}/assets/empty_mito_depth.tsv", checkIfExists: true)

    ch_mitogenome_assembly_results = buildAssemblyUploadRows(
        ch_canonical_upload_rows,
        ch_variant_upload_rows,
        ch_mito_depth.map { meta, tsv -> [ meta.mt_assembly_prefix, tsv ] },
        no_depth_file,
        params.mitogenome_summary_min_length
    )

    //
    // SUBWORKFLOW: UPLOAD_RESULTS
    //
    // Conditional uploading of results to SQL and species check - only run if not skipped
    // All these processes access the OceanOmics PostgreSQL database.
    def ch_qc_input = Channel.empty()
    if (!params.skip_upload_results && params.sql_config) {
        // Per-sample circularity-check evidence from both assemblers feeds the QC
        // gate (anomaly block; the circular condition itself comes via meta.circular).
        // Keyed by mt_assembly_prefix and total over the canonical assemblies, which is what
        // lets the gate attach it with a plain join instead of collecting the whole channel.
        ch_mitogenome_circularity_evidence = ch_canonical_circularity_evidence

        UPLOAD_RESULTS (
            ch_mitogenome_assembly_results,
            ch_mitogenome_annotation_results,
            ch_mitogenome_blast_results,
            ch_mitogenome_lca_results,
            ch_mitogenome_lca_raw_results,
            ch_mitogenome_circularity_evidence,
            ch_mitogenome_region_counts,
            sql_config // params.sql_config
        )

        // Keyed on mt_assembly_prefix, not the whole meta map. Both sides descend from
        // MITOGENOME_ANNOTATION so their metas usually agree, but each is restored from its own
        // cache entry on a -resume, and a whole-map join that stops agreeing drops every sample
        // SILENTLY -- the failure mode that emptied the QC gate one operator upstream (see the
        // upload-receipt join in upload_results_mito). Hardened defensively: mt_assembly_prefix
        // is a field of both metas, so prefix equality is implied by map equality and this can
        // only ever match the same pairs or more, and it is 1:1 on both sides (one QC row and
        // one annotation bundle per assembly) so there is no fan-out.
        ch_qc_input = UPLOAD_RESULTS.out.qc_ready
            .map { meta, species_name, proceed_qc, circular ->
                [ meta.mt_assembly_prefix, meta, species_name, proceed_qc, circular ]
            }
            .join(
                ch_mitogenome_annotation_results.map { meta, files -> [ meta.mt_assembly_prefix, files ] },
                by: 0
            )
            .map { _prefix, meta, species_name, proceed_qc, circular, annotation_files ->
                [ meta, species_name, proceed_qc, circular, annotation_files ]
            }

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
