/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../../nf-core/utils_nfcore_pipeline'

// Mitogenome assembly
include { GETORGANELLE_CONFIG       } from '../../../../modules/nf-core/getorganelle/config'
include { CAT_FASTQ                 } from '../../../../modules/nf-core/cat/fastq'
include { FASTP                     } from '../../../../modules/nf-core/fastp'
include { GETORGANELLE_FROMREADS    } from '../../../../modules/nf-core/getorganelle/fromreads'
include { GETORGANELLE_RESEED       } from '../../../../modules/local/getorganelle/reseed'
include { GETORGANELLE_GENEDB       } from '../../../../modules/local/getorganelle/genedb'
include { MITOHIFI_FINDMITOREFERENCE } from '../../../../modules/nf-core/mitohifi/findmitoreference'
include { RELABEL_REFERENCE_GB      } from '../../../../modules/local/relabel_reference_gb'
include { GETORGANELLE_JOIN         } from '../../../../modules/local/getorganelle/join'
include { GETORGANELLE_CHECK        } from '../../../../modules/local/getorganelle/check'
include { PUSH_MTDNA_ASSM_RESULTS   } from '../../../../modules/local/upload_results/mtdna'

// Read the corrected circular verdict from a GETORGANELLE_CHECK evidence TSV.
// Returns true / false (final_verdict_circular) or null when the column is
// missing / NA / unparseable. Folded back into meta.circular so the non-circular
// scaffolds the reference test confirms circular are annotated (MITOS2) and gated
// (GenBank QC) as circular. Defined at file scope so it resolves inside .map closures.
def parseFinalVerdictCircular(tsv) {
    try {
        def rows = tsv.text.readLines()
        if (rows.size() < 2) return null
        def header = rows[0].split('\t')
        def idx = header.findIndexOf { it.trim() == 'final_verdict_circular' }
        if (idx < 0) return null
        def cells = rows[1].split('\t')
        if (idx >= cells.size()) return null
        def v = cells[idx].trim().toLowerCase()
        return (v == 'true') ? true : (v == 'false' ? false : null)
    } catch (ignored) {
        return null
    }
}

// Decide whether a first-pass GetOrganelle result warrants a reseed attempt.
// A circular genome is the only outright success: anything else (empty FASTA,
// a fragmented multi-contig result, or a single non-circular contig) is worth
// retrying with a closely-related seed BEFORE we fall back to concatenating
// scaffolds just to scrape a species ID. The reseed result is never re-assessed
// here, so this gate fires at most once per sample.
def needsReseed(fasta, logFile) {
    // No assembled contig at all -> always worth a reseed.
    try { if (fasta == null || fasta.size() == 0) return true } catch (ignored) { return true }
    // More than one contig -> fragmented assembly, worth a reseed.
    try {
        if (fasta.text.readLines().count { it.startsWith('>') } > 1) return true
    } catch (ignored) { /* fall through to the log-based check below */ }
    def txt
    try { txt = logFile.text } catch (ignored) { return false }
    // A circular genome is the success case; any non-circular result reseeds.
    if (txt =~ /Result status of .*: circular genome/) return false
    return true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ASSEMBLY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ASSEMBLY_GETORG {

    take:
    fastp_reads // tuple val(meta), path(fastp)
    organelle_type // "animal_mt" passed to download the correct congig files
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_summary_files = Channel.empty()
    // Per-variant assembly inputs for the DB upload. Each GetOrganelle variant
    // (first-pass, reseed, reference-guided-join) is captured as its own
    // mitogenome_data row, keyed by its variant prefix (the fasta basename).
    // Entries are [ variant_prefix, meta, fasta, log ].
    ch_db_variant_inputs = Channel.empty()
    // Per-sample reference GenBank (findMitoReference, via RELABEL_REFERENCE_GB).
    // Only reseed candidates download one, so this is a partial channel; the
    // annotation subworkflow falls back to a fresh lookup / curated asset for the
    // rest. Stays empty when the reseed stage is skipped.
    ch_reference_gb = Channel.empty()

     
    //
    // MODULE: Set up for running GetOrganelle
    //

    GETORGANELLE_CONFIG (
        organelle_type // val(organelle_type)
    )

    // Split fastp_reads based on whether concatenation is needed
    fastp_reads_split = fastp_reads.branch { meta, reads ->
        def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
        def needsConcatenation = meta.single_end ? readList.size > 1 : readList.size > 2
        
        needs_concat: needsConcatenation
            return [meta, reads]
        no_concat: !needsConcatenation
            return [meta, reads]
    }
    //
    // MODULE: Concatenate fastq reads where there are multiple fastq R1 and R2 files
    //

    CAT_FASTQ (
        fastp_reads_split.needs_concat
    )
    
    // Combine the results
    final_reads = fastp_reads_split.no_concat.mix(CAT_FASTQ.out.reads)

    //
    // MODULE: Trim raw HiC reads (adapter + poly-G) before assembly.
    //   HiC reads enter this pipeline raw -- the draft-genome pipeline only trims
    //   its Illumina short reads -- so residual Illumina adapter + NovaSeq poly-G
    //   tails get assembled onto contig ends and block circularisation. fastp
    //   strips them. Illumina (ilmn) reads were already fastp-trimmed upstream, so
    //   they bypass this step to avoid double-trimming. Gated by skip_hic_fastp.
    //
    reads_for_assembly = final_reads
    if (!params.skip_hic_fastp) {
        reads_by_type = final_reads.branch { meta, _reads ->
            hic:   meta.sequencing_type == 'hic'
            other: true
        }
        FASTP (
            reads_by_type.hic,
            []   // adapter_fasta: empty -> fastp auto-detects adapters by overlap
        )
        reads_for_assembly = reads_by_type.other.mix(FASTP.out.reads)

        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { _meta, json -> json })
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.tool_params.collect { it[1] })
        ch_versions      = ch_versions.mix(FASTP.out.versions.first())
    }

    //
    // Get the GetOrganelle version to be used in the naming. Ensures the proper version is always used.
    //

    version_ch = GETORGANELLE_CONFIG.out.versions
    .map { versions_file ->
        def content = versions_file.text
        def pattern = /getorganelle:\s*([^\s\n\r]+)/
        def matcher = content =~ pattern
        return matcher ? matcher[0][1].replaceAll(/["']/, '') : "unknown"
    }

    //
    // Embed the assembly prefix into meta.
    //

    fastp_with_mt_assembly_prefix = reads_for_assembly.combine(version_ch)
    .map { meta, fasta, version ->  // Destructure all 3 elements correctly
        def version_stripped = version.replaceAll('\\.', '')
        def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.getorg${version_stripped}"
        def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
        [meta_ext, fasta]  // Return only what you want
    }


    //
    // Combine fastp files with the GetOrganelle database
    //   .first() makes the single-item db channel a reusable value channel so it
    //   can be combined into both the first pass and the reseed pass.
    //

    ch_db = GETORGANELLE_CONFIG.out.db.first()
    combined_input = fastp_with_mt_assembly_prefix.combine(ch_db)

    //
    // MODULE: Run assembly using GetOrganelle from reads (first pass, generic seed)
    //

    GETORGANELLE_FROMREADS (
        combined_input // tuple val(meta), path(fastp), val(organelle_type), path(db)
    )

    //
    // PHASE 3: Conditional reseed.
    //   Assess each first-pass result; for poor/failed assemblies (seed too
    //   divergent), download a closely related mitogenome and re-run seeded.
    //   Gated by params.skip_getorganelle_reseed (default: run).
    //

    // Join first-pass fasta + log + reads so we can assess and, if needed, reseed.
    ch_firstpass = GETORGANELLE_FROMREADS.out.fasta
        .join(GETORGANELLE_FROMREADS.out.log, by: 0)
        .join(fastp_with_mt_assembly_prefix, by: 0)   // [meta, fasta, log, reads]

    if (!params.skip_getorganelle_reseed) {

        ch_assessed = ch_firstpass.branch { _meta, fasta, log, _reads ->
            reseed: needsReseed(fasta, log)
            keep:   !needsReseed(fasta, log)
        }

        // Split reseed candidates by taxon. INVERTEBRATES (corals) seed from the
        // curated coral mitogenome DB -- a broad, label-free Anthozoa seed + label
        // database, so a wrong/coarse species label can no longer pick a wrong seed.
        // VERTEBRATES keep the per-sample findMitoReference download. The coral DB is
        // NEVER used for vertebrates.
        ch_reseed_branched = ch_assessed.reseed.branch { meta, _fasta, _log, _reads ->
            invert: meta.invertebrates
            vert:   true
        }

        // --- Vertebrate reseed path: findMitoReference seed + custom gene db. ---
        MITOHIFI_FINDMITOREFERENCE (
            ch_reseed_branched.vert.map { meta, _fasta, _log, _reads -> meta }
        )
        ch_vert_seed = MITOHIFI_FINDMITOREFERENCE.out.reference
            .filter { _meta, ref_fasta, ref_gb -> ref_fasta.size() > 0 && ref_gb.size() > 0 }
            .map { meta, ref_fasta, _ref_gb -> [meta, ref_fasta] }   // [meta, seed]
        ch_vert_ref_gb = MITOHIFI_FINDMITOREFERENCE.out.reference
            .filter { _meta, ref_fasta, ref_gb -> ref_fasta.size() > 0 && ref_gb.size() > 0 }
            .map { meta, _ref_fasta, ref_gb -> [meta, ref_gb] }      // [meta, ref_gb]

        // Custom label database from the reference GenBank (disentangles contigs for
        // divergent animal mitogenomes); the vertebrate reseed only runs WITH it.
        GETORGANELLE_GENEDB ( ch_vert_ref_gb )

        // Relabel the reference GenBank to a per-sample name for the assembly summary.
        RELABEL_REFERENCE_GB ( ch_vert_ref_gb )
        // Per-sample reference GenBank (vertebrates only). Coral annotation now picks
        // its reference from the DB by sequence, so corals need none here.
        ch_reference_gb = RELABEL_REFERENCE_GB.out.gb

        // --- Invertebrate (coral) reseed path: broad coral DB seed + label db. ---
        ch_coral_seed  = Channel.fromPath("${projectDir}/assets/coral_mito_refdb.fasta",       checkIfExists: true).first()
        ch_coral_label = Channel.fromPath("${projectDir}/assets/coral_mito_refdb.label.fasta", checkIfExists: true).first()
        ch_invert_seed  = ch_reseed_branched.invert
            .map { meta, _fasta, _log, _reads -> meta }.combine(ch_coral_seed)
            .map { meta, seed -> [meta, seed] }
        ch_invert_genes = ch_reseed_branched.invert
            .map { meta, _fasta, _log, _reads -> meta }.combine(ch_coral_label)
            .map { meta, genes -> [meta, genes] }

        // Merge the two reseed paths. Inverts always carry both seed + genes (assets);
        // verts carry them only when findMitoReference + GENEDB succeeded.
        ch_seed  = ch_vert_seed.mix(ch_invert_seed)
        ch_genes = GETORGANELLE_GENEDB.out.genes.mix(ch_invert_genes)

        ch_reseed_reads = ch_assessed.reseed.map { meta, _fasta, _log, reads -> [meta, reads] }

        // Reseed candidates that got BOTH a seed AND a usable gene database. The inner
        // joins drop any sample missing either input, enforcing "never reseed without
        // --genes" (only ever bites a vertebrate with a sparse/absent reference).
        ch_reseed_input = ch_reseed_reads
            .join(ch_seed, by: 0)                       // [meta, reads, seed]
            .join(ch_genes, by: 0)                      // [meta, reads, seed, genes]
            .combine(ch_db)                             // [meta, reads, seed, genes, org_type, db]
            .map { meta, reads, seed, genes, org_type, db -> [meta, reads, org_type, db, seed, genes] }

        GETORGANELLE_RESEED (
            ch_reseed_input // tuple val(meta), path(fastp), val(organelle_type), path(db), path(seed), path(genes)
        )

        // Reseed candidates that could NOT be reseeded (no reference seed, or a
        // reference too sparsely annotated to build a gene database): fall back to
        // their first-pass result so no sample is dropped. Keyed off "ready"
        // samples (those with both a seed and a gene database).
        ch_reseed_ready = ch_seed.join(ch_genes, by: 0)
            .map { meta, _seed, _genes -> [meta, true] }   // [meta, true]

        ch_reseed_fallback = ch_assessed.reseed
            .map { meta, fasta, log, _reads -> [meta, fasta, log] }
            .join(ch_reseed_ready, by: 0, remainder: true)
            .filter { it[3] == null }                      // not ready -> fall back
            .map { meta, fasta, log, _ready -> [meta, fasta, log] }

        ch_reseed_seedless     = ch_reseed_fallback.map { meta, fasta, _log -> [meta, fasta] }
        ch_reseed_seedless_log = ch_reseed_fallback.map { meta, _fasta, log -> [meta, log] }

        // Resolve each seeded reseed candidate to the better of (reseed result,
        // first-pass result). Prefer the reseed output, but if it came back empty
        // keep the first-pass assembly so the multi-contig concat / species-ID
        // fallback (SANITISE_FASTA) still has something to work with.
        ch_reseed_firstpass = ch_assessed.reseed
            .map { meta, fasta, log, _reads -> [meta, fasta, log] }   // [meta, fp_fasta, fp_log]

        ch_reseed_resolved = GETORGANELLE_RESEED.out.fasta
            .join(GETORGANELLE_RESEED.out.log, by: 0)                 // [meta, rs_fasta, rs_log]
            .join(ch_reseed_firstpass, by: 0)                        // [meta, rs_fasta, rs_log, fp_fasta, fp_log]
            .map { meta, rs_fasta, rs_log, fp_fasta, fp_log ->
                def useReseed = rs_fasta && rs_fasta.size() > 0
                [meta, useReseed ? rs_fasta : fp_fasta, useReseed ? rs_log : fp_log]
            }

        // Final per-sample assembly = kept first-pass + resolved reseed + seedless fallback.
        ch_assembly_fasta = ch_assessed.keep.map { meta, fasta, _log, _reads -> [meta, fasta] }
            .mix(ch_reseed_resolved.map { meta, fasta, _log -> [meta, fasta] })
            .mix(ch_reseed_seedless)
        ch_assembly_log = ch_assessed.keep.map { meta, _fasta, log, _reads -> [meta, log] }
            .mix(ch_reseed_resolved.map { meta, _fasta, log -> [meta, log] })
            .mix(ch_reseed_seedless_log)

        // Reseed assembly as its own DB variant (keyed by fasta basename).
        ch_db_variant_inputs = ch_db_variant_inputs.mix(
            GETORGANELLE_RESEED.out.fasta
                .join(GETORGANELLE_RESEED.out.log, by: 0)
                .map { meta, fasta, log -> [ fasta.baseName, meta, fasta, log ] }
        )

        // Reseed-specific files to fold into the collected channels below.
        ch_multiqc_files = ch_multiqc_files
            .mix(GETORGANELLE_RESEED.out.tool_params.collect { it[1] })
            .mix(GETORGANELLE_GENEDB.out.tool_params.collect { it[1] })
            .mix(MITOHIFI_FINDMITOREFERENCE.out.tool_params.collect { it[1] })
        ch_summary_files = ch_summary_files
            .mix(GETORGANELLE_RESEED.out.fasta.map { _meta, fasta -> fasta })
            .mix(GETORGANELLE_RESEED.out.log.map { _meta, log -> log })
            .mix(GETORGANELLE_RESEED.out.org_assm_graph)
            .mix(GETORGANELLE_RESEED.out.raw_assm_graph)
            .mix(GETORGANELLE_RESEED.out.simp_assm_graph)
            .mix(RELABEL_REFERENCE_GB.out.gb.map { _meta, gb -> gb })
        ch_versions = ch_versions
            .mix(GETORGANELLE_RESEED.out.versions.first())
            .mix(GETORGANELLE_GENEDB.out.versions.first())
            .mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())

    } else {
        ch_assembly_fasta = GETORGANELLE_FROMREADS.out.fasta
        ch_assembly_log   = GETORGANELLE_FROMREADS.out.log
    }

    //
    // Collect files
    //

    // ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_FROMREADS.out.etc.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_CONFIG.out.tool_params)
    ch_multiqc_files = ch_multiqc_files.mix(CAT_FASTQ.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_FROMREADS.out.tool_params.collect { it[1] })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.fasta.map { meta, fasta -> fasta })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.log.map { meta, log -> log })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.org_assm_graph)
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.raw_assm_graph)
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_FROMREADS.out.simp_assm_graph)
    // First-pass assembly as its own DB variant (keyed by fasta basename).
    ch_db_variant_inputs = ch_db_variant_inputs.mix(
        GETORGANELLE_FROMREADS.out.fasta
            .join(GETORGANELLE_FROMREADS.out.log, by: 0)
            .map { meta, fasta, log -> [ fasta.baseName, meta, fasta, log ] }
    )
    ch_versions = ch_versions.mix(GETORGANELLE_CONFIG.out.versions.first())
    ch_versions = ch_versions.mix(GETORGANELLE_FROMREADS.out.versions.first())
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // Capture assembly topology (circular vs linear) into meta so downstream
    // annotation (MITOS2) can choose the right mode: GetOrganelle assembles a
    // circular molecule linearised at an arbitrary point, so running MITOS with
    // --linear severs whichever gene straddles that break. GetOrganelle reports
    // "Result status of ...: circular genome" in its log (the same signal used
    // by needsReseed above); absence means non-circular, a failed parse means
    // unknown (null). Enrich fasta + log together so both keep an identical meta
    // map — the parent workflow joins these two channels by the whole meta map.
    //
    ch_assembly_topo = ch_assembly_fasta
        .join(ch_assembly_log, by: 0)
        .map { meta, fasta, log ->
            def circ
            try { circ = (log && (log.text =~ /Result status of .*: circular genome/)) ? true : false }
            catch (ignored) { circ = null }
            [ meta + [ circular: circ ], fasta, log ]
        }
    ch_assembly_fasta = ch_assembly_topo.map { meta, fasta, _log -> [ meta, fasta ] }
    ch_assembly_log   = ch_assembly_topo.map { meta, _fasta, log -> [ meta, log ] }

    //
    // MODULE: Re-test circularity for non-circular GetOrganelle scaffolds + screen
    //   length/tandem-repeat anomalies. A single scaffold GetOrganelle could not
    //   formally close is frequently a complete circle linearised elsewhere; the
    //   check confirms this against the related reference (full reference coverage)
    //   and corrects the circular verdict. Samples with no findMitoReference get an
    //   empty placeholder reference (the check then just records "no_reference").
    //
    // Join the reference by a stable key (mt_assembly_prefix): ch_reference_gb's
    // meta predates the meta.circular enrichment above, so a whole-meta join would
    // never match. ch_reference_gb is partial (reseed/vertebrate samples only), so
    // remainder + placeholder keeps every fasta sample.
    def no_reference_gb = file("${projectDir}/assets/NO_REFERENCE.gb", checkIfExists: true)
    // Pair every assembly with its reference (or null), keyed by mt_assembly_prefix
    // (ch_reference_gb's meta predates the meta.circular enrichment, so a whole-meta
    // join would never match; it is also partial, so remainder keeps every fasta).
    // ch_reference_gb is single-use, so this is its only consumer; the reference is
    // then carried forward (through GETORGANELLE_JOIN's passthrough) rather than
    // re-joined.
    ch_ref_keyed = ch_reference_gb.map { meta, ref -> [ meta.mt_assembly_prefix, ref ] }
    ch_fasta_ref = ch_assembly_fasta
        .map { meta, fasta -> [ meta.mt_assembly_prefix, meta, fasta ] }
        .join(ch_ref_keyed, by: 0, remainder: true)
        .filter { items -> items[1] != null }   // keep fasta rows; drop reference-only remainder
        .map { items ->
            def ref = items.size() > 3 ? items[3] : null
            [ items[1], items[2], ref ]          // [meta, fasta, ref-or-null]
        }

    //
    // MODULE: Reference-guided scaffold join. A multi-scaffold GetOrganelle result
    //   that still carries the full gene set (a control-region repeat the reseed
    //   could not span) is ordered/oriented against the related reference and joined
    //   into one molecule (named "<prefix>_rgj"), so GETORGANELLE_CHECK and the
    //   annotation see a single record instead of a blind concatenation. Only
    //   multi-scaffold samples WITH a real reference are joined; single-scaffold
    //   results and multi-scaffold results without a reference pass through (the
    //   latter are concatenated downstream by SANITISE_FASTA -> "<prefix>_concat").
    //
    ch_join_branched = ch_fasta_ref.branch { _meta, fasta, ref ->
        def multi = false
        try { multi = fasta.text.readLines().count { it.startsWith('>') } > 1 } catch (ignored) { multi = false }
        rgj:    multi && ref != null
        direct: true
    }

    GETORGANELLE_JOIN (
        ch_join_branched.rgj
    )

    // Re-pair the joined assembly with its (passed-through) reference, mix back the
    // pass-through samples (substituting the empty placeholder where no reference
    // exists), and feed the combined channel to the circularity/length check.
    ch_getorg_check_in = GETORGANELLE_JOIN.out.fasta
        .join(GETORGANELLE_JOIN.out.reference, by: 0)
        .map { meta, joined, ref -> [ meta, joined, ref ] }
        .mix( ch_join_branched.direct.map { meta, fasta, ref -> [ meta, fasta, ref ?: no_reference_gb ] } )

    // The joined assembly replaces the multi-scaffold fasta for every downstream
    // stage (annotation, LCA, QC, upload), not just the check.
    ch_assembly_fasta = ch_getorg_check_in.map { meta, fasta, _ref -> [ meta, fasta ] }

    GETORGANELLE_CHECK (
        ch_getorg_check_in
    )

    // Join evidence + the joined "_rgj" fasta feed the summary / multiqc alongside
    // the check evidence. The fasta is what gives the "_rgj" assembly its
    // final_length_bp / contig count in the assembly summary (the pre-join
    // scaffolds carry a different prefix); without it the rgj row has no length.
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_JOIN.out.fasta.map { _meta, fasta -> fasta })
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_JOIN.out.evidence.map { _meta, ev -> ev })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_JOIN.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(GETORGANELLE_JOIN.out.versions.first())

    // RGJ (reference-guided-join) assembly as its own DB variant. It produces no
    // GetOrganelle log of its own, so reuse the pre-join resolved log (same reads
    // -> same coverage); ch_assembly_log here is still the pre-check version.
    ch_db_variant_inputs = ch_db_variant_inputs.mix(
        GETORGANELLE_JOIN.out.fasta
            .join(ch_assembly_log, by: 0)
            .map { meta, fasta, log -> [ fasta.baseName, meta, fasta, log ] }
    )

    // Fold the corrected circular verdict back into meta on BOTH the fasta and log
    // channels (the parent workflow joins them by the whole meta map, so the two
    // must stay identical). null verdict -> leave meta unchanged.
    ch_circular_update = GETORGANELLE_CHECK.out.evidence
        .map { meta, ev -> [ meta, parseFinalVerdictCircular(ev) ] }

    ch_assembly_fasta = ch_assembly_fasta
        .join(ch_circular_update, by: 0)
        .map { meta, fasta, v -> [ (v == null ? meta : meta + [ circular: v ]), fasta ] }
    ch_assembly_log = ch_assembly_log
        .join(ch_circular_update, by: 0)
        .map { meta, logf, v -> [ (v == null ? meta : meta + [ circular: v ]), logf ] }

    // Per-variant DB upload results. The corrected circular verdict only applies to
    // the variant that was actually checked (the final assembly); its evidence file
    // is named after that variant's prefix, so key the verdict by prefix and attach
    // it to the matching variant only. Every other variant carries no verdict, so
    // the DB push falls back to that variant's own GetOrganelle log topology.
    ch_checked_circ = GETORGANELLE_CHECK.out.evidence
        .map { _meta, ev -> [ ev.name.replaceAll(/\.getorg_check\.tsv$/, ''), parseFinalVerdictCircular(ev) ] }

    ch_db_assembly_results = ch_db_variant_inputs
        .map { prefix, meta, fasta, logf -> [ prefix, [ meta, fasta, logf ] ] }
        .join( ch_checked_circ, by: 0, remainder: true )
        .filter { it[1] != null }                       // drop check-only remainders
        .map { items ->
            def prefix  = items[0]
            def payload = items[1]
            def verdict = items.size() > 2 ? items[2] : null
            def meta    = payload[0] + [ mt_assembly_prefix: prefix ]
            if (verdict != null) meta = meta + [ circular: verdict ]
            [ meta, payload[1], payload[2] ]
        }

    // Evidence feeds the assembly summary (anomaly reason + circular override) and
    // is emitted for the QC gate (anomaly block + circular condition).
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_CHECK.out.evidence.map { _meta, ev -> ev })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_CHECK.out.evidence.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_CHECK.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(GETORGANELLE_CHECK.out.versions.first())

    //
    // Emit outputs
    //

    emit:
    assembly_fasta  = ch_assembly_fasta
    assembly_log    = ch_assembly_log
    db_assembly_results = ch_db_assembly_results   // channel: [ meta(per-variant prefix+circular), fasta, log ] -> one DB row per variant
    reference_gb    = ch_reference_gb              // channel: [ meta, reference.gb ] (partial: reseed candidates only)
    circularity_evidence = GETORGANELLE_CHECK.out.evidence  // channel: [ meta, getorg_check.tsv ]
    summary_files   = ch_summary_files
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
