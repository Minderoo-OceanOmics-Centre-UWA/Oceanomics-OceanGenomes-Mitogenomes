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
include { REFERENCE_DIVERGENCE      } from '../../../../modules/local/reference_divergence'
include { REFERENCE_CANDIDATES      } from '../../../../modules/local/reference_candidates'
include { REFERENCE_RANK            } from '../../../../modules/local/reference_rank'
include { GETORGANELLE_JOIN         } from '../../../../modules/local/getorganelle/join'
include { GETORGANELLE_CHECK        } from '../../../../modules/local/getorganelle/check'

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

// Reduce the per-variant DB inputs to the PROVENANCE variants -- every assembly attempt this
// sample made EXCEPT the one that won, which the parent workflow already emits as the sample's
// canonical row. Excluding the winner here is what lets the parent plain-mix the two channels
// instead of reconciling them with a whole-run groupTuple. See the call site for the run that
// cost.
//
// The exclusion used to be "everything except the FIRST PASS", on the reasoning that the
// first-pass basename equalled meta.mt_assembly_prefix and so collided with the canonical row.
// That held only while the canonical row was misnamed: it was filed under the sample-level
// prefix even when the molecule in it was a reseed. Now the canonical row carries the winning
// assembly's own name, so the row that would collide is the WINNER's, and the first pass is
// exactly the provenance worth keeping -- a superseded attempt with its own name, its own
// stats and no depth.
//
// variant_inputs is [ variant_prefix, meta, fasta, log ]; checked_circ is
// [ run_prefix, checked_identity, verdict ], one row per sample.
//
// combine(by: 0) on the LINEAGE key, not a join(remainder: true) on identity: every sample
// runs exactly one GETORGANELLE_CHECK, so the key always matches and each variant emits as
// soon as its OWN sample is checked. A remainder join would hold every unmatched variant --
// which is all of them but the checked one -- until the channel closed.
def selectProvenanceVariants(variant_inputs, checked_circ) {
    return variant_inputs
        .map { prefix, meta, fasta, logf -> [ meta.mt_assembly_run_prefix, prefix, meta, fasta, logf ] }
        .combine( checked_circ, by: 0 )
        .filter { _run_prefix, prefix, _meta, _fasta, _logf, checked_identity, _verdict ->
            prefix != checked_identity
        }
        .map { _run_prefix, prefix, meta, fasta, logf, _checked_identity, _verdict ->
            // A superseded attempt never takes the corrected circular verdict: that verdict
            // describes the molecule that was checked, which by definition is not this one.
            // Each falls back to its own GetOrganelle log topology in the DB push.
            [ meta + [ mt_assembly_prefix: prefix ], fasta, logf ]
        }
}

// Read the tier from a REFERENCE_DIVERGENCE flag file. Mirrors the helper in the
// MitoHiFi subworkflow; both routes grade the resolved reference the same way.
def parseReferenceTier(tsv) {
    try {
        return tsv.text.readLines().find { it?.trim() }?.split('\t', -1)?.first()?.trim()?.toUpperCase() ?: 'UNKNOWN'
    } catch (ignored) {
        return 'UNKNOWN'
    }
}

// Did REFERENCE_RANK actually choose a reference? Mirrors the helper in the MitoHiFi
// subworkflow. It emits EMPTY chosen_reference files when it declined to substitute
// (no candidate recruited any reads, no reads subsampled, or a lone unmapped
// candidate), so an evidence-free substitution never displaces the reference
// findMitoReference resolved. OG2021 and OG810 are the samples this catches on the
// GetOrganelle side. Truthiness alone is not enough; the files exist either way.
def hasChosenReference(fasta, gb) {
    return fasta && gb && fasta.size() > 0 && gb.size() > 0
}

// Did REFERENCE_CANDIDATES return any candidate references for this sample? Read from its
// always-emitted status TSV (status column == 'found'), so re-selection can branch on a
// per-sample value instead of on the presence of the optional candidates output -- which
// could only be classified once the channel closed. Defined at file scope so it resolves
// inside the .map closures.
def referenceCandidatesFound(statusFile) {
    try {
        def rows = statusFile.text.readLines().findAll { it?.trim() }
        if (rows.size() < 2) return false
        def idx = rows[0].split('\t').findIndexOf { it.trim() == 'status' }
        if (idx < 0) return false
        def cells = rows[1].split('\t', -1)
        return idx < cells.size() && cells[idx].trim() == 'found'
    } catch (ignored) {
        return false
    }
}

// Did GETORGANELLE_GENEDB build a usable gene database for this sample? Read from its
// always-emitted status TSV (ready column == 'yes'), so the reseed fallback can branch on a
// per-sample value rather than a remainder join against the optional genes output. A sample
// whose reference was too sparsely annotated reports ready=no and falls back to its first pass.
def genedbReady(statusFile) {
    try {
        def rows = statusFile.text.readLines().findAll { it?.trim() }
        if (rows.size() < 2) return false
        def idx = rows[0].split('\t').findIndexOf { it.trim() == 'ready' }
        if (idx < 0) return false
        def cells = rows[1].split('\t', -1)
        return idx < cells.size() && cells[idx].trim().toLowerCase() == 'yes'
    } catch (ignored) {
        return false
    }
}

// Is this divergence tier worth re-selecting a reference for? CONGENERIC already has
// the best obtainable reference and UNKNOWN carries no evidence the reference is
// poor. CROSS_ORDER is not special-cased here: GetOrganelle seeds from the reference
// rather than recruiting reads against it, so it degrades gracefully where MitoHiFi
// fails outright, and re-selection is the right response for it too.
def shouldReselectReference(tier) {
    return tier in ['CONFAMILIAL', 'DIFFERENT_FAMILY', 'NON_CONGENERIC', 'CROSS_ORDER']
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
    // rest. Stays empty when the reseed stage is skipped. This is the EMITTED
    // channel and stays partial on purpose -- REFERENCE_RELEVANCE (downstream) only
    // grades samples that actually resolved a reference.
    ch_reference_gb = Channel.empty()
    // Explicit NO_REFERENCE.gb placeholder rows ([meta, no_reference_gb]) for every
    // sample that will NEVER resolve a real reference (kept first pass, coral reseed,
    // seedless vertebrate, or the whole run when reseed is skipped). Emitted at the
    // moment the branch decision is made so the GETORGANELLE_CHECK reference join below
    // can be a plain per-sample join instead of a remainder join that waits on the
    // whole run to close. Kept OUT of the emitted reference_gb above.
    ch_reference_placeholders = Channel.empty()
    // The empty placeholder reference. Hoisted here (from the check-join site) so both
    // reseed branches AND the reseed-skipped else-branch can mix it in.
    def no_reference_gb = file("${projectDir}/assets/placeholders/NO_REFERENCE.gb", checkIfExists: true)

     
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
    // Embed the assembly prefix into meta, as TWO fields with different lifetimes.
    //
    //   mt_assembly_prefix     -- IDENTITY. What this molecule is called. Re-stamped to the
    //                             FASTA basename every time curation renames the assembly
    //                             (reseed, _rgj, _collapsed, _concat), so it always names the
    //                             molecule in hand. Drives publishing, DB rows and every
    //                             post-annotation join.
    //   mt_assembly_run_prefix -- LINEAGE. This sample's GetOrganelle run. Written ONCE, here,
    //                             and never reassigned. Drives the joins that reunite artefacts
    //                             which are per sample+assembler and identical across variants
    //                             (the reads, the resolved reference, this sample's mtdna bundle).
    //
    // They were one field until a curated variant proved they are not the same thing: the
    // assembly stage kept the sample-level value while annotation overwrote it with the FASTA
    // basename, so every join spanning that rewrite silently missed. That dropped 23 of 168
    // assemblies at the QC gate and filed the reseed's depth against the first-pass row.
    // Anything reassigning mt_assembly_run_prefix is a bug; the identity is the one that moves.
    //
    // Deliberately NOT set on anything upstream of here (GETORGANELLE_CONFIG, CAT_FASTQ,
    // FASTP): read prep predates the version lookup this prefix needs, and keying it on
    // meta.id / the samplesheet assembly_prefix keeps those tasks out of the hash change.

    fastp_with_mt_assembly_prefix = reads_for_assembly.combine(version_ch)
    .map { meta, fasta, version ->  // Destructure all 3 elements correctly
        def version_stripped = version.replaceAll('\\.', '')
        def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.getorg${version_stripped}"
        def meta_ext = meta + [
            mt_assembly_prefix:     mt_assembly_prefix,
            mt_assembly_run_prefix: mt_assembly_prefix
        ]
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

        // Split reseed candidates by taxon. INVERTEBRATES seed from the curated
        // per-phylum mitogenome DB for their class (assets/refdb/<group>/) -- a broad,
        // label-free seed + label database, so a wrong/coarse species label can no
        // longer pick a wrong seed. VERTEBRATES keep the per-sample findMitoReference
        // download. The invert DBs are NEVER used for vertebrates.
        ch_reseed_branched = ch_assessed.reseed.branch { meta, _fasta, _log, _reads ->
            invert: meta.invertebrates
            vert:   true
        }

        // --- Vertebrate reseed path: findMitoReference seed + custom gene db. ---
        MITOHIFI_FINDMITOREFERENCE (
            ch_reseed_branched.vert.map { meta, _fasta, _log, _reads -> meta }
        )
        ch_vert_reference = MITOHIFI_FINDMITOREFERENCE.out.reference
            .filter { _meta, ref_fasta, ref_gb -> ref_fasta.size() > 0 && ref_gb.size() > 0 }

        //
        // Reference divergence guard + re-selection for the vertebrate reseed.
        // The reseed seeds GetOrganelle from this reference AND builds the custom
        // gene label DB from it, so a reference that fell back past genus degrades
        // both. Grade it, and when it is not congeneric replace it with the candidate
        // the sample's own reads map to best. Until now the GetOrganelle route had no
        // divergence flag at all, which is why its assemblies showed an empty tier in
        // the assembly summary.
        //
        REFERENCE_DIVERGENCE (
            ch_vert_reference.map { meta, _ref_fasta, ref_gb -> [meta, ref_gb] }
        )

        ch_vert_reads = ch_reseed_branched.vert.map { meta, _fasta, _log, reads -> [meta, reads] }

        ch_vert_route = ch_vert_reference
            .join(ch_vert_reads, by: 0)                       // [meta, ref_fasta, ref_gb, reads]
            .join(REFERENCE_DIVERGENCE.out.flag, by: 0)       // + flag
            .branch { _meta, _ref_fasta, _ref_gb, _reads, flag ->
                reselect: params.enable_reference_reselection &&
                          shouldReselectReference(parseReferenceTier(flag))
                keep: true
            }

        REFERENCE_CANDIDATES (
            ch_vert_route.reselect.map { meta, _ref_fasta, _ref_gb, _reads, _flag -> meta }
        )

        // Branch the reselect samples on whether REFERENCE_CANDIDATES actually returned any
        // candidates -- read from its always-emitted status (a per-sample VALUE), not from the
        // presence/absence of the optional candidates output. This turns the old remainder join
        // (which could not classify a no-candidate sample until the whole run closed) into a
        // plain per-sample branch, so re-selection no longer holds the reseed fork open.
        ch_reselect_branched = ch_vert_route.reselect
            .join(REFERENCE_CANDIDATES.out.status.map { meta, s -> [ meta, referenceCandidatesFound(s) ] }, by: 0)
            .branch { _meta, _ref_fasta, _ref_gb, _reads, _flag, found ->
                rank: found
                keep: true
            }

        REFERENCE_RANK (
            ch_reselect_branched.rank
                .map { meta, _ref_fasta, _ref_gb, reads, _flag, _found -> [meta, reads] }
                .join(REFERENCE_CANDIDATES.out.candidates, by: 0)
        )

        // Resolve each vertebrate reseed to its reference. REFERENCE_RANK.out.reference is total
        // over the `rank` branch (RANK runs on exactly those samples), so a plain join suffices;
        // an empty chosen_reference means RANK declined to substitute, so keep the original.
        // Samples with no candidates keep their original reference, and the keep-path reference
        // is untouched. No remainder, no channel-close dependency.
        ch_vert_resolved = ch_vert_route.keep
            .mix( ch_reselect_branched.keep
                    .map { meta, ref_fasta, ref_gb, reads, flag, _found -> [meta, ref_fasta, ref_gb, reads, flag] } )
            .mix( ch_reselect_branched.rank
                    .map { meta, ref_fasta, ref_gb, reads, flag, _found -> [meta, ref_fasta, ref_gb, reads, flag] }
                    .join(REFERENCE_RANK.out.reference, by: 0)
                    .map { meta, ref_fasta, ref_gb, reads, flag, chosen_fasta, chosen_gb ->
                        hasChosenReference(chosen_fasta, chosen_gb)
                            ? [ meta, chosen_fasta, chosen_gb, reads, flag ]
                            : [ meta, ref_fasta, ref_gb, reads, flag ]
                    } )

        ch_vert_seed = ch_vert_resolved
            .map { meta, ref_fasta, _ref_gb, _reads, _flag -> [meta, ref_fasta] }   // [meta, seed]
        ch_vert_ref_gb = ch_vert_resolved
            .map { meta, _ref_fasta, ref_gb, _reads, _flag -> [meta, ref_gb] }      // [meta, ref_gb]

        // Custom label database from the reference GenBank (disentangles contigs for
        // divergent animal mitogenomes); the vertebrate reseed only runs WITH it.
        GETORGANELLE_GENEDB ( ch_vert_ref_gb )

        // Relabel the reference GenBank to a per-sample name for the assembly summary.
        RELABEL_REFERENCE_GB ( ch_vert_ref_gb )
        // Per-sample reference GenBank (vertebrates only). Coral annotation now picks
        // its reference from the DB by sequence, so corals need none here.
        ch_reference_gb = RELABEL_REFERENCE_GB.out.gb

        // Totality for the GETORGANELLE_CHECK reference join: emit an explicit
        // NO_REFERENCE.gb placeholder (keyed by lineage like the real references) for
        // every sample that will never carry a RELABEL_REFERENCE_GB row --
        //   * kept first passes (never reseed),
        //   * coral reseeds (seed the curated coral DB; never RELABEL),
        //   * vertebrate reseeds whose findMitoReference found nothing (no ref at all).
        // Vertebrates WITH a reference but no gene database still carry the real RELABEL
        // row (they fall back to first-pass via ch_reseed_seedless yet keep their
        // resolved reference), so they are deliberately NOT placeheld here -- doing so
        // would double-key them and duplicate the sample at the plain join. The
        // vert-no-ref set is the exact complement of the non-empty filter at :311.
        ch_reference_placeholders = ch_assessed.keep
            .map { meta, _fasta, _log, _reads -> [ meta, no_reference_gb ] }
            .mix( ch_reseed_branched.invert.map { meta, _fasta, _log, _reads -> [ meta, no_reference_gb ] } )
            .mix( MITOHIFI_FINDMITOREFERENCE.out.reference
                    .filter { _meta, ref_fasta, ref_gb -> !(ref_fasta.size() > 0 && ref_gb.size() > 0) }
                    .map { meta, _ref_fasta, _ref_gb -> [ meta, no_reference_gb ] } )

        // --- Invertebrate reseed path: the curated DB for the sample's phylum. ---
        // InvertTaxonGroups.seedDbGroup() resolves the class to a group directory under
        // assets/refdb/, or to null when no curated database covers that class. A null
        // is NOT backfilled with some other phylum's database: every invertebrate used
        // to be reseeded from the Anthozoa DB, so a mollusc or a sea star was re-run
        // against a seed far too divergent to assemble from, failing a second time after
        // burning the full GetOrganelle walltime. An unmapped class is simply not
        // reseeded and keeps its first-pass assembly (via ch_reseed_readiness below).
        ch_invert_group = ch_reseed_branched.invert
            .map { meta, _fasta, _log, _reads ->
                def group = InvertTaxonGroups.seedDbGroup(meta.class)
                if (group) {
                    log.info "GETORGANELLE_RESEED: ${meta.id} (class ${meta.class}) seeding from assets/refdb/${group}"
                } else {
                    log.warn "GETORGANELLE_RESEED: ${meta.id} has no curated seed database for class " +
                             "'${meta.class}' -- keeping the first-pass assembly instead of reseeding " +
                             "from another phylum. Add the class to InvertTaxonGroups.seedDbGroup() and " +
                             "build its database with bin/build_invert_reference_db.py to reseed it."
                }
                [ meta, group ]
            }

        ch_invert_seeded = ch_invert_group.filter { _meta, group -> group != null }
        ch_invert_seed  = ch_invert_seeded.map { meta, group ->
            [ meta, file("${projectDir}/assets/refdb/${group}/${group}_mito_refdb.fasta", checkIfExists: true) ] }
        ch_invert_genes = ch_invert_seeded.map { meta, group ->
            [ meta, file("${projectDir}/assets/refdb/${group}/${group}_mito_refdb.label.fasta", checkIfExists: true) ] }

        // Merge the two reseed paths. Inverts carry both seed + genes (assets) whenever
        // their class resolves to a curated database; verts carry them only when
        // findMitoReference + GENEDB succeeded.
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

        // Reseed candidates that could NOT be reseeded (no reference seed, or a reference too
        // sparsely annotated to build a gene database): fall back to their first-pass result so
        // no sample is dropped. Readiness is now a per-sample VALUE, total over the whole reseed
        // set and disjoint across its three sources, so the fallback is a plain join + filter
        // rather than a remainder join that could not classify a not-ready sample until the run
        // closed:
        //   * inverts: ready iff their class resolved to a curated seed database -> its group
        //   * vertebrates with a reference: ready iff GETORGANELLE_GENEDB built  -> its status
        //   * vertebrates with no findMitoReference at all                       -> not ready
        // (a vertebrate with a reference but a sparse gene database reports ready=no from GENEDB,
        // so it falls back here while still keeping its resolved reference upstream.)
        ch_reseed_readiness = ch_invert_group
                .map { meta, group -> [ meta, group != null ] }
            .mix( GETORGANELLE_GENEDB.out.status.map { meta, s -> [ meta, genedbReady(s) ] } )
            .mix( MITOHIFI_FINDMITOREFERENCE.out.reference
                    .filter { _meta, ref_fasta, ref_gb -> !(ref_fasta.size() > 0 && ref_gb.size() > 0) }
                    .map { meta, _ref_fasta, _ref_gb -> [ meta, false ] } )

        ch_reseed_fallback = ch_assessed.reseed
            .map { meta, fasta, log, _reads -> [meta, fasta, log] }
            .join(ch_reseed_readiness, by: 0)
            .filter { _meta, _fasta, _log, ready -> !ready }   // not ready -> fall back
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
            // Divergence tier + re-selection audit trail, same as the MitoHiFi route.
            .mix(REFERENCE_DIVERGENCE.out.flag.map { _meta, flag -> flag })
            .mix(REFERENCE_RANK.out.ranking.map { _meta, ranking -> ranking })
            .mix(REFERENCE_CANDIDATES.out.status.map { _meta, status -> status })
        ch_versions = ch_versions
            .mix(GETORGANELLE_RESEED.out.versions.first())
            .mix(GETORGANELLE_GENEDB.out.versions.first())
            .mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())
            .mix(REFERENCE_DIVERGENCE.out.versions.first())
            .mix(REFERENCE_CANDIDATES.out.versions.first())
            .mix(REFERENCE_RANK.out.versions.first())

    } else {
        ch_assembly_fasta = GETORGANELLE_FROMREADS.out.fasta
        ch_assembly_log   = GETORGANELLE_FROMREADS.out.log
        // No reseed => no reference is ever resolved; give every sample an explicit
        // placeholder so the check-join below stays a plain per-sample join instead of
        // a remainder join that waits for the whole run to close.
        ch_reference_placeholders = GETORGANELLE_FROMREADS.out.fasta
            .map { meta, _fasta -> [ meta, no_reference_gb ] }
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
    // This is also where the assembly's IDENTITY is settled for the reseed fork. By here
    // ch_assembly_fasta holds whichever molecule won -- the kept first pass (basename
    // "<run>") or the reseed result (basename "<run>reseed") -- so stamping
    // mt_assembly_prefix from the FASTA basename is what makes the curated variant name
    // itself. Doing it here rather than at the three mixes above is deliberate: the log's
    // basename is not the assembly's name, so the identity can only be derived once fasta
    // and log are joined and can be stamped onto both at once.
    ch_assembly_topo = ch_assembly_fasta
        .join(ch_assembly_log, by: 0)
        .map { meta, fasta, log ->
            def circ
            try { circ = (log && (log.text =~ /Result status of .*: circular genome/)) ? true : false }
            catch (ignored) { circ = null }
            [ meta + [ circular: circ, mt_assembly_prefix: fasta.baseName ], fasta, log ]
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
    // Join the reference by the LINEAGE key: one reference is resolved per sample and it is
    // the same reference whichever variant won, so it is not per-assembly. It also has to be
    // a key that survives curation -- the reference meta predates both the meta.circular
    // enrichment and the identity re-stamp above, so a whole-meta join would never match and
    // an identity join would miss every renamed assembly.
    //
    // Real references (RELABEL_REFERENCE_GB, partial) mixed with the per-sample placeholders
    // emitted at each branch decision above make ch_ref_keyed TOTAL over ch_assembly_fasta --
    // exactly one row per sample, no duplicates -- so a plain join classifies every sample the
    // moment its own reference/placeholder arrives, instead of waiting for both channels to
    // close (the remainder-join whole-run barrier this replaces). The placeholder set is built
    // from channels already local at each branch, so establishing totality costs no join and
    // introduces no channel-close dependency. Single-use: the reference is then carried forward
    // (through GETORGANELLE_JOIN's passthrough) rather than re-joined.
    ch_ref_keyed = ch_reference_gb.mix(ch_reference_placeholders)
        .map { meta, ref -> [ meta.mt_assembly_run_prefix, ref ] }
    ch_fasta_ref = ch_assembly_fasta
        .map { meta, fasta -> [ meta.mt_assembly_run_prefix, meta, fasta ] }
        .join(ch_ref_keyed, by: 0)
        .map { _key, meta, fasta, ref -> [ meta, fasta, ref ] }   // [meta, fasta, ref]

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
        // ref is now always a file (the NO_REFERENCE.gb placeholder when none resolved),
        // so "has a real reference" is ref.size() > 0, not ref != null.
        rgj:    multi && ref.size() > 0
        direct: true
    }

    GETORGANELLE_JOIN (
        ch_join_branched.rgj
    )

    // Re-pair the joined assembly with its (passed-through) reference, mix back the
    // pass-through samples (substituting the empty placeholder where no reference
    // exists), and feed the combined channel to the circularity/length check.
    //
    // Re-stamp identity here too: GETORGANELLE_JOIN renames the molecule to "<prefix>_rgj",
    // and this channel feeds BOTH the downstream assembly and GETORGANELLE_CHECK, so the
    // check's evidence has to carry the joined assembly's name -- that identity is what
    // decides, downstream, which variant the corrected circular verdict belongs to and which
    // DB row the QC gate attaches the evidence to. Pass-through samples restamp to the same
    // value they already had.
    ch_getorg_check_in = GETORGANELLE_JOIN.out.fasta
        .join(GETORGANELLE_JOIN.out.reference, by: 0)
        .map { meta, joined, ref -> [ meta, joined, ref ] }
        .mix( ch_join_branched.direct.map { meta, fasta, ref -> [ meta, fasta, ref ?: no_reference_gb ] } )
        .map { meta, fasta, ref -> [ meta + [ mt_assembly_prefix: fasta.baseName ], fasta, ref ] }

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

    // Carry the "_rgj" identity re-stamp onto the log as well, now that the DB-variant mix
    // above (which needs the PRE-stamp meta, since GETORGANELLE_JOIN.out.fasta carries it) is
    // done. An rgj assembly produces no GetOrganelle log of its own, so the log's own basename
    // is not the assembly's name and cannot be stamped from it -- re-attach it to the stamped
    // meta by LINEAGE, which is stable across the rename. Without this the fasta and log stop
    // sharing a meta map and the whole-meta joins just below, and in the parent workflow,
    // silently drop every _rgj sample.
    ch_assembly_log = ch_assembly_log
        .map { meta, logf -> [ meta.mt_assembly_run_prefix, logf ] }
        .join(ch_assembly_fasta.map { meta, _fasta -> [ meta.mt_assembly_run_prefix, meta ] }, by: 0)
        .map { _key, logf, meta -> [ meta, logf ] }

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

    // Per-variant DB upload results, for the PROVENANCE variants only -- see
    // selectProvenanceVariants above for which variant is excluded and why.
    //
    // This tuple carries BOTH kinds of key, because it does two different jobs:
    //   [0] LINEAGE  -- fans this one per-sample check out to all of that sample's variants.
    //   [1] IDENTITY -- names the assembly the verdict actually describes, so only that one
    //                   takes it. A superseded first pass and the reseed that replaced it
    //                   share a lineage but are different molecules with different topology.
    //
    // The identity used to be recovered by parsing the evidence FILENAME, because meta had
    // nowhere to carry it -- mt_assembly_prefix was the sample-level value at this point. It
    // is stamped from the checked FASTA's basename now, so read it from meta and drop the
    // parse: one less place where a name has to be taken apart to recover something the
    // pipeline already knew.
    ch_checked_circ = GETORGANELLE_CHECK.out.evidence
        .map { meta, ev ->
            [ meta.mt_assembly_run_prefix,
              meta.mt_assembly_prefix,
              parseFinalVerdictCircular(ev) ]
        }

    ch_db_assembly_results = selectProvenanceVariants(ch_db_variant_inputs, ch_checked_circ)

    // Evidence feeds the assembly summary (anomaly reason + circular override) and
    // is emitted for the QC gate (anomaly block + circular condition).
    ch_summary_files = ch_summary_files.mix(GETORGANELLE_CHECK.out.evidence.map { _meta, ev -> ev })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_CHECK.out.evidence.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_CHECK.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(GETORGANELLE_CHECK.out.versions.first())

    // Per-sample bundle of everything this stage publishes into <prefix>/mtdna
    // (final assembly + its GetOrganelle log + the circularity/anomaly check), keyed by the
    // LINEAGE prefix -- it is this sample's assembly-run folder, not any one molecule's, and
    // the key has to survive the collapse rename so the mirror can restage the whole folder
    // into <prefix>_collapsed/mtdna for genuinely collapsed samples.
    // Sized with groupKey so each sample's bundle releases as soon as its own three files
    // (final assembly + GetOrganelle log + circularity/anomaly check) arrive, rather than
    // waiting for the whole run to close. The count is EXACTLY 3 for every sample: all three
    // channels derive from ch_getorg_check_in (one row each per sample), so no group is ever
    // short and a plain groupTuple releases every sample -- no close-time remainder is needed
    // or wanted. Unwrap the GroupKey back to the plain prefix so the collapse mirror's key
    // join still matches.
    ch_mtdna_files = ch_assembly_fasta
        .mix( ch_assembly_log,
              GETORGANELLE_CHECK.out.evidence )
        .map { meta, f -> [ groupKey(meta.mt_assembly_run_prefix, 3), f ] }
        .groupTuple()
        .map { key, files -> [ key.getGroupTarget(), files ] }

    //
    // Emit outputs
    //

    emit:
    mtdna_files     = ch_mtdna_files               // channel: [ mt_assembly_run_prefix, [ mtdna files ] ]
    assembly_fasta  = ch_assembly_fasta
    assembly_log    = ch_assembly_log
    // CONTRACT: the PROVENANCE variants only -- every attempt this sample made EXCEPT the one
    // that won, whose identity is the parent workflow's canonical row for that assembly. The
    // parent relies on that to mix this channel straight into the upload candidates without a
    // de-duplicating groupTuple (which would be a whole-run barrier). The exclusion is by
    // IDENTITY (see selectProvenanceVariants), so it holds however curation renamed the winner.
    db_assembly_results = ch_db_assembly_results   // channel: [ meta(per-variant identity), fasta, log ] -> one DB row per superseded attempt
    reference_gb    = ch_reference_gb              // channel: [ meta, reference.gb ] (partial: reseed candidates only)
    // CONTRACT: exactly one row per assembly FASTA this subworkflow emits. That holds by
    // construction here -- assembly_fasta and GETORGANELLE_CHECK are both derived from
    // ch_getorg_check_in, so every emitted assembly is checked and no placeholder is needed
    // (unlike the MitoHiFi subworkflow, which has failed / no-contig branches to fill in).
    // The parent workflow and the QC gate rely on that totality to attach evidence with a
    // plain per-key join rather than collecting the whole channel; if a future change lets an
    // assembly bypass the check, it must emit assets/placeholders/empty_circularity_check.tsv for that
    // sample or the sample will be dropped at the QC gate.
    circularity_evidence = GETORGANELLE_CHECK.out.evidence  // channel: [ meta, getorg_check.tsv ]
    // Reads keyed by the LINEAGE prefix, for MITOGENOME_COVERAGE in the parent workflow.
    // Keyed rather than meta-joined because meta gains `circular` and a re-stamped identity
    // downstream, so a whole-meta join would never match. Lineage rather than identity because
    // these are READS, not an assembly: there is exactly one set per sample+assembler and every
    // GetOrganelle variant is built from precisely those, so one remap covers whichever variant
    // wins and the reads have no curated name to be keyed by. Which MOLECULE gets measured is
    // decided on the other side of that join -- the parent passes the sanitised assembly, i.e.
    // the curated winner -- and the resulting depth is filed under that molecule's identity.
    depth_reads     = fastp_with_mt_assembly_prefix.map { m, r -> [ m.mt_assembly_run_prefix, r ] }
    summary_files   = ch_summary_files
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
