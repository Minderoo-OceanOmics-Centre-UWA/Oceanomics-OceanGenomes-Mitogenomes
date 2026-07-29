/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../../nf-core/utils_nfcore_pipeline'

// Mitogenome assembly
include { MITOHIFI_FINDMITOREFERENCE       } from '../../../../modules/nf-core/mitohifi/findmitoreference'
include { CAT_FASTQ                         } from '../../../../modules/nf-core/cat/fastq'
include { MITOHIFI_MITOHIFI    } from '../../../../modules/nf-core/mitohifi/mitohifi'
include { MITOHIFI_AVERAGE_COVERAGE        } from '../../../../modules/local/mitohifi/average_coverage'
include { MITOHIFI_CHECK_CIRCULARITY       } from '../../../../modules/local/mitohifi/check_circularity'
include { RELABEL_REFERENCE_GB             } from '../../../../modules/local/relabel_reference_gb'
include { REFERENCE_DIVERGENCE             } from '../../../../modules/local/reference_divergence'
include { REFERENCE_CANDIDATES             } from '../../../../modules/local/reference_candidates'
include { REFERENCE_RANK                   } from '../../../../modules/local/reference_rank'
include { OATK                             } from '../../../../modules/local/oatk'
include { OATK_CHECK                       } from '../../../../modules/local/oatk/check_circularity'
include { ASSEMBLY_NO_RESULT               } from '../../../../modules/local/assembly_no_result'
include { PUSH_MTDNA_ASSM_RESULTS   } from '../../../../modules/local/upload_results/mtdna'

// Read the circularity verdict from a MITOHIFI_CHECK_CIRCULARITY evidence TSV.
// Returns true / false (final_verdict_circular) or null when the column is
// missing / NA / unparseable. Folded into meta.circular so the annotation stage
// (EMMA / MITOS2) and the GenBank QC gate see the real HiFi topology, mirroring
// the GetOrganelle path. Defined at file scope so it resolves inside .map closures.
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

// Is this divergence tier worth re-selecting a reference for? CONGENERIC is not (a
// same-genus reference is already the best obtainable), CROSS_ORDER is routed to
// reference-free assembly instead, and UNKNOWN carries no evidence that the
// reference is poor -- so only the explicitly non-congeneric tiers qualify.
def shouldReselectReference(tier) {
    return tier in ['CONFAMILIAL', 'DIFFERENT_FAMILY', 'NON_CONGENERIC']
}

// Count CDS features in MitoHiFi's own final_mitogenome.gb. MitoHiFi annotates the
// assembly it produces, so the protein-coding-gene count is available here, inside
// the assembly subworkflow -- no dependency on the downstream annotation
// subworkflow, and therefore no dataflow cycle when the count is used to route.
// Returns null when the file is missing/unparseable, which routing treats as "no
// evidence of a collapse" so an odd GenBank never diverts a good assembly.
def countGenbankCds(gb) {
    try {
        if (!gb || !gb.exists() || gb.size() == 0) return null
        return gb.text.readLines().count { it.startsWith('     CDS ') }
    } catch (ignored) {
        return null
    }
}

def parseReferenceTier(tsv) {
    try {
        return tsv.text.readLines().find { it?.trim() }?.split('\t', -1)?.first()?.trim()?.toUpperCase() ?: 'UNKNOWN'
    } catch (ignored) {
        return 'UNKNOWN'
    }
}

def parseReferenceStatus(tsv) {
    try {
        def rows = tsv.text.readLines()
        if (rows.size() < 2) return 'lookup_error'
        def header = rows[0].split('\t', -1)
        def idx = header.findIndexOf { it.trim() == 'status' }
        def values = rows[1].split('\t', -1)
        return idx >= 0 && idx < values.size() ? values[idx].trim().toLowerCase() : 'lookup_error'
    } catch (ignored) {
        return 'lookup_error'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ASSEMBLY WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ASSEMBLY_MITOHIFI {

    take:
    fastp_reads // tuple val(meta), path(fastp)
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_summary_files = Channel.empty()

    def mitohifi_version = params.mitohifi_container.tokenize(':').last()
    def mitohifi_version_stripped = mitohifi_version.replaceAll('\\.', '')

    //
    // map just the meta for the species query
    //

    ch_species_reference = fastp_reads
    .map { meta, _files ->
        def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${mitohifi_version_stripped}mitohifi"
        meta + [ mt_assembly_prefix: mt_assembly_prefix ]
    }
    // .view()
    
    //
    // MODULE: Find a closely related species for reference
    //

    MITOHIFI_FINDMITOREFERENCE (
        ch_species_reference
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
    // MODULE: Concatenate fastq reads where there are multiple fastq files
    //

    CAT_FASTQ (
        fastp_reads_split.needs_concat
    )

    // Combine the results
    final_reads = fastp_reads_split.no_concat.mix(CAT_FASTQ.out.reads)

    //
    // Combine fastp files with the mito reference output
    //

    // Embed the assembly prefix before the reference join so reads remain routable
    // when findMitoReference emits a no-reference placeholder. The version is parsed from the
    // pinned container tag (params.mitohifi_container) rather than `mitohifi.py
    // --version`, because the 3.2.3 release ships a stale self-reported version
    // (3.2.1). The tag is the single source of truth: bump it in nextflow.config and
    // the version in every assembly name follows automatically.
    //

    ch_reads_by_prefix = final_reads.map { meta, reads ->
        def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${mitohifi_version_stripped}mitohifi"
        def meta_ext = meta + [ mt_assembly_prefix: mt_assembly_prefix ]
        [meta_ext, reads]
    }

    ch_reference_outcomes = MITOHIFI_FINDMITOREFERENCE.out.reference
        .join(MITOHIFI_FINDMITOREFERENCE.out.status, by: 0)
    ch_reference_joined = ch_reads_by_prefix.join(ch_reference_outcomes, by: 0)
    ch_reference_branched = ch_reference_joined.branch { _meta, _reads, ref_fasta, ref_gb, _status ->
        found: ref_fasta.size() > 0 && ref_gb.size() > 0
        missing: true
    }
    combined_with_mt_assembly_prefix = ch_reference_branched.found
        .map { meta, reads, ref_fasta, ref_gb, _status -> [meta, reads, ref_fasta, ref_gb] }

    //
    // MODULE: Pre-assembly reference divergence guard.
    // findMitoReference walks the sample's NCBI lineage and grabs the first
    // complete mitogenome, so a species with no congeneric record (deep-sea / poorly
    // sampled taxa) silently gets a divergent reference. MitoHiFi then drops the most
    // divergent gene blocks during reference-based read recruitment, yielding a
    // clean-looking but gene-incomplete collapse. Compare sample vs reference
    // taxonomy up front and record a CONGENERIC/.../NON_CONGENERIC review flag; the
    // assembly summary turns anything non-congeneric into a manual_review reason.
    // Taxonomy-only, always exits 0.
    //
    // Runs on the raw reference GenBank rather than the relabelled one: relabelling
    // now happens *after* re-selection, so that the name (and everything downstream
    // that reads species/accession from it) describes the reference actually used.
    //

    REFERENCE_DIVERGENCE (
        combined_with_mt_assembly_prefix.map { meta, _reads, _ref_fasta, ref_gb -> [meta, ref_gb] }
    )

    //
    // Route on the divergence tier.
    //   CROSS_ORDER      -> skip MitoHiFi entirely; recruitment against a reference
    //                       from another order maps almost nothing. Reference-free.
    //   not congeneric   -> re-select: fetch several candidates and keep the one the
    //                       sample's own reads map to best.
    //   congeneric       -> keep it. A same-genus reference is the best obtainable,
    //                       so there is nothing to re-select and no call to spend.
    // UNKNOWN (unparseable taxonomy) is deliberately left on the keep path: without a
    // reliable tier there is no evidence the reference is poor, and re-selecting on a
    // guess would churn NCBI for every sample with a thin lineage.
    //
    ch_reference_route = combined_with_mt_assembly_prefix
        .join(REFERENCE_DIVERGENCE.out.flag, by: 0)
        .branch { _meta, _reads, _ref_fasta, _ref_gb, flag ->
            oatk_direct: parseReferenceTier(flag) == 'CROSS_ORDER'
            reselect: params.enable_reference_reselection &&
                      shouldReselectReference(parseReferenceTier(flag))
            keep: true
        }

    //
    // MODULE: Reference re-selection (REFERENCE_CANDIDATES -> REFERENCE_RANK).
    // Replaces findMitoReference's first-hit reference with the candidate that
    // recruits the most of this sample's reads. See modules/local/reference_rank.
    //
    REFERENCE_CANDIDATES (
        ch_reference_route.reselect.map { meta, _reads, _ref_fasta, _ref_gb, _flag -> meta }
    )

    REFERENCE_RANK (
        ch_reference_route.reselect
            .map { meta, reads, _ref_fasta, _ref_gb, _flag -> [meta, reads] }
            .join(REFERENCE_CANDIDATES.out.candidates, by: 0)
    )

    // Substitute the chosen reference. remainder:true keeps samples that produced no
    // candidates (NCBI lookup failed, or nothing usable came back) on their original
    // reference, so re-selection can only ever improve on the previous behaviour.
    ch_reselected = ch_reference_route.reselect
        .join(REFERENCE_RANK.out.reference, by: 0, remainder: true)
        .filter { it[1] != null }   // drop any right-only remainder
        .map { items ->
            def chosen_fasta = items.size() > 5 ? items[5] : null
            def chosen_gb    = items.size() > 6 ? items[6] : null
            (chosen_fasta && chosen_gb)
                ? [ items[0], items[1], chosen_fasta, chosen_gb, items[4] ]
                : [ items[0], items[1], items[2], items[3], items[4] ]
        }

    ch_reference_resolved = ch_reference_route.keep.mix(ch_reselected)

    //
    // MODULE: Relabel the reference GenBank to a per-sample name so it can be fed
    // (collision-free) to the assembly summary, which reads the reference species
    // and accession from it. Runs on the resolved reference, so the summary and the
    // relevance check report the reference the assembly was actually built from.
    //

    RELABEL_REFERENCE_GB (
        ch_reference_resolved.map { meta, _reads, _ref_fasta, ref_gb, _flag -> [meta, ref_gb] }
    )


    //
    // MODULE: Run assembly using MitoHifi from reads
    //

    MITOHIFI_MITOHIFI (
        ch_reference_resolved.map { meta, reads, ref_fasta, ref_gb, _flag -> [meta, reads, ref_fasta, ref_gb] },
        "r",
        "2"  // Fallback genetic code only: the module prefers meta.genetic_code (derived per-sample from taxonomic class in PREPARE_SAMPLESHEET). 2 = vertebrate mitochondrial.
    )

    // Branch the assembled fasta on emptiness: when MitoHiFi finishes without
    // producing a final mitogenome the wrapper emits an empty placeholder.
    // Only run the coverage step for samples that actually assembled.
    ch_mitohifi_fasta_branched = MITOHIFI_MITOHIFI.out.fasta.branch { meta, fasta ->
        assembled: fasta.size() > 0
        failed:    fasta.size() == 0
    }

    ch_average_coverage_input = ch_mitohifi_fasta_branched.assembled
        .join(MITOHIFI_MITOHIFI.out.stats, by: 0)
        .join(MITOHIFI_MITOHIFI.out.coverage_mapping, by: 0)
        .map { meta, _fasta, stats, cov_map -> [meta, stats, cov_map] }

    MITOHIFI_AVERAGE_COVERAGE (
        ch_average_coverage_input
    )

    //
    // MODULE: Re-test circularity for assemblies MitoHiFi flagged non-circular.
    //   MitoHiFi's terminal-overlap check yields false negatives on hifiasm
    //   assemblies that are genuinely circular (the closed unitig loses its
    //   self-overlap once MitoHiFi rotates/trims it). Left uncorrected, the false
    //   flag is written to the SQL db as a "scaffold" and trips the assembly
    //   summary's not_circularised manual-review reason. The check remaps the
    //   already-mapped HiFi reads to a doubled reference (junction-spanning reads)
    //   and reads the hifiasm c/l contig flag, then corrects was_circular in the
    //   with-coverage stats. The corrected table is a drop-in replacement (same
    //   basename) so the SQL upload and summary pick it up with no further changes.
    //

    // gb is joined for the control-region location of any length-inflating tandem
    // repeat. MitoHiFi writes final_mitogenome.gb whenever it writes the final
    // FASTA, so every assembled sample carries one and the inner join drops none.
    ch_check_circularity_input = ch_mitohifi_fasta_branched.assembled
        .join(MITOHIFI_AVERAGE_COVERAGE.out.stats, by: 0)
        .join(MITOHIFI_MITOHIFI.out.coverage_mapping, by: 0)
        .join(MITOHIFI_MITOHIFI.out.gb, by: 0)
        .map { meta, fasta, stats, cov_map, gb -> [meta, stats, fasta, cov_map, gb] }

    MITOHIFI_CHECK_CIRCULARITY (
        ch_check_circularity_input
    )

    // Corrected with-coverage stats supersede the average-coverage output for
    // every downstream consumer (SQL upload, assembly summary, MultiQC).
    ch_assembled_stats = MITOHIFI_CHECK_CIRCULARITY.out.stats

    // For samples that failed to assemble, fall back to the (empty)
    // contigs_stats.tsv as a log placeholder so PUSH_MTDNA_ASSM_RESULTS still
    // gets a tuple to process.
    ch_failed_assembly_log = ch_mitohifi_fasta_branched.failed
        .join(MITOHIFI_MITOHIFI.out.stats, by: 0)
        .map { meta, _fasta, stats -> [meta, stats] }

    ch_assembly_log = ch_assembled_stats.mix(ch_failed_assembly_log)

    //
    // MODULE: OATK reference-free fallback (gated by params.enable_oatk_fallback).
    //   MitoHiFi recruits reads by mapping to a related-species reference; for taxa
    //   with no close NCBI relative (a divergent / cross-order reference) that maps
    //   ~zero reads and hifiasm produces no contig. Oatk instead identifies the
    //   mitogenome by profile-HMM over a de-novo HiFi assembly, so it needs no
    //   species reference and recovers exactly these divergent-reference failures.
    //   Runs on the three ways a reference can fail a sample: no reference resolved
    //   at all, a cross-order reference (routed before assembly), and a reference
    //   distant enough that MitoHiFi returned a gene-incomplete collapse rather than
    //   nothing. Emits a FASTA that the parent workflow feeds into annotation. Off by
    //   default (needs an Oatk container + an OatkDB mito profile for the sample clade).
    //
    ch_oatk_fasta = Channel.empty()
    ch_oatk_log = Channel.empty()
    ch_oatk_circularity_evidence = Channel.empty()
    ch_routed_failure_fasta = Channel.empty()
    ch_routed_failure_log = Channel.empty()
    // Oatk's reads keyed by the OATK assembly prefix. Kept separate from
    // ch_reads_by_prefix because ch_oatk_input below OVERWRITES mt_assembly_prefix
    // with the oatk prefix, so a prefix-keyed join against ch_reads_by_prefix would
    // match nothing and silently drop every oatk sample from the depth measurement.
    ch_oatk_reads = Channel.empty()

    ch_direct_oatk_reads = ch_reference_branched.missing
        .map { meta, reads, _ref_fasta, _ref_gb, status ->
            def outcome = parseReferenceStatus(status)
            def reason = outcome == 'lookup_error' ? 'reference_lookup_error' : 'no_reference'
            [meta, reads, reason]
        }
        .mix(ch_reference_route.oatk_direct
            .map { meta, reads, _ref_fasta, _ref_gb, _flag -> [meta, reads, 'cross_order_reference'] })
    if (params.enable_oatk_fallback) {
        // Reads keyed by the SAME enriched meta (with mt_assembly_prefix) the failed
        // branch carries, so the join matches. combined_with_mt_assembly_prefix holds
        // [meta_ext, reads, ref_fasta, ref_gb]; take meta + reads.
        // Oatk assembler tag for the run prefix, parsed from the container version
        // (1.0 -> v10oatk) so the summary files it as an Oatk run distinct from the
        // failed MitoHiFi attempt. Overwrite mt_assembly_prefix so every downstream
        // output (Oatk contig, annotation, summary) is named and grouped under it.
        def oatk_version_stripped = (params.oatk_container ?: 'oatk:1.0')
            .tokenize(':').last().tokenize('--').first().replaceAll('\\.', '')

        ch_failed_oatk_reads = ch_mitohifi_fasta_branched.failed
            .join(ch_reads_by_prefix, by: 0)
            .map { meta, _empty_fasta, reads -> [meta, reads, 'empty_mitohifi'] }

        // An empty FASTA is not the only way a divergent reference ruins an assembly.
        // When the reference is too distant, MitoHiFi's reference-guided read
        // recruitment drops the most divergent gene blocks and returns a clean-looking
        // but gene-incomplete molecule: OG2102 (Rouleina attrita, non-congeneric
        // Alepocephalus reference) came back as a plausible 15.6 kb contig carrying
        // 7 of 13 protein-coding genes, so it never reached the empty-FASTA fallback.
        // Route those to Oatk too -- it needs no reference, so it is exactly the tool
        // for a collapse the reference caused. The PCG count comes from MitoHiFi's own
        // final_mitogenome.gb (inner join: an assembly with no GenBank yields no
        // evidence and is left alone).
        def expected_pcg_count = (params.mitogenome_summary_expected_pcg_count ?: 13) as int
        ch_gene_incomplete_oatk_reads = ch_mitohifi_fasta_branched.assembled
            .join(MITOHIFI_MITOHIFI.out.gb, by: 0)
            .map { meta, _fasta, gb -> [meta, countGenbankCds(gb)] }
            .filter { _meta, cds -> cds != null && cds < expected_pcg_count }
            .join(ch_reads_by_prefix, by: 0)
            .map { meta, _cds, reads -> [meta, reads, 'gene_incomplete_mitohifi'] }

        // The MitoHiFi assembly is deliberately NOT withdrawn when this fires: it stays
        // published alongside the Oatk attempt (distinct v10oatk prefix), so the summary
        // shows both and curation picks. Discarding it would lose information.
        ch_oatk_input = ch_direct_oatk_reads
            .mix(ch_failed_oatk_reads, ch_gene_incomplete_oatk_reads)
            .map { meta, reads, reason ->
                def oatk_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${oatk_version_stripped}oatk"
                [ meta + [ mt_assembly_prefix: oatk_prefix, circular: null, assembler_fallback: 'oatk', fallback_reason: reason ], reads ]
            }

        // Same reads, now carrying the oatk prefix, for the uniform depth measurement.
        ch_oatk_reads = ch_oatk_input

        // Stage the .fam AND its .h3* nhmmer indexes together (nhmmscan needs them
        // side by side); params.oatk_mito_db points at the .fam.
        ch_oatk_db = Channel.fromPath("${params.oatk_mito_db}*", checkIfExists: true).collect()

        OATK (
            ch_oatk_input,
            ch_oatk_db
        )

        // Run the same FASTA-vs-reference circularity/anomaly check the other assemblers
        // get, on assembled oatk contigs only (skip the no-contig empties). Without it the
        // oatk contig carries no circularity evidence, so in the parent workflow it can
        // only leave the collapse-concatemer join via the (channel-close) remainder path,
        // i.e. it reaches annotation as a delayed end-of-run batch instead of incrementally
        // like MitoHiFi/GetOrganelle. Emitting matching evidence (identical meta) makes it a
        // matched join item that flows straight through, and gives it a real meta.circular
        // for MITOS2 and the GenBank QC gate.
        ch_oatk_assembled = OATK.out.fasta.filter { _meta, fasta -> fasta.size() > 0 }

        // Assembled oatk contigs always emit a GFA; pair them for the self-link
        // circularity read.
        ch_oatk_fa_gfa = ch_oatk_assembled.join(OATK.out.gfa, by: 0)   // [meta, fasta, gfa]

        // Attach the per-sample reference GenBank. RELABEL_REFERENCE_GB.out.gb carries the
        // MitoHiFi-prefixed meta (not the oatk one), so key on [id, sequencing_type, date].
        // The true no-reference oatk path (ch_reference_branched.missing) has no relabelled
        // reference, so remainder + the NO_REFERENCE placeholder keeps every oatk contig
        // (check_getorganelle.py treats a zero-length reference as absent).
        def no_reference_gb_oatk = file("${projectDir}/assets/NO_REFERENCE.gb", checkIfExists: true)
        ch_oatk_ref_keyed = RELABEL_REFERENCE_GB.out.gb
            .map { m, gb -> [ [m.id, m.sequencing_type, m.date], gb ] }
        ch_oatk_check_in = ch_oatk_fa_gfa
            .map { m, fasta, gfa -> [ [m.id, m.sequencing_type, m.date], m, fasta, gfa ] }
            .join(ch_oatk_ref_keyed, by: 0, remainder: true)
            .filter { it[1] != null }   // keep oatk rows; drop reference-only remainder
            .map { items ->
                def ref = (items.size() > 4 && items[4] != null) ? items[4] : no_reference_gb_oatk
                [ items[1], items[2], items[3], ref ]   // [meta, fasta, gfa, ref]
            }

        OATK_CHECK ( ch_oatk_check_in )

        // Fold the corrected circular verdict into meta on every oatk channel that is later
        // joined by the whole meta map (fasta <-> evidence in the parent collapse join; the
        // QC-gate evidence). All three must carry the SAME meta, so fold the identical
        // verdict into each. remainder:true keeps the no-contig oatk (no evidence ->
        // circular:null, empty FASTA filtered out of annotation downstream).
        ch_oatk_circ_verdict = OATK_CHECK.out.evidence
            .map { m, tsv -> [ m, parseFinalVerdictCircular(tsv) ] }

        ch_oatk_fasta = OATK.out.fasta
            .join(ch_oatk_circ_verdict, by: 0, remainder: true)
            .map { m, fasta, circ -> [ m + [ circular: circ ], fasta ] }
        ch_oatk_log = OATK.out.log
            .join(ch_oatk_circ_verdict, by: 0, remainder: true)
            .map { m, log, circ -> [ m + [ circular: circ ], log ] }
        ch_oatk_circularity_evidence = OATK_CHECK.out.evidence
            .map { m, ev -> [ m + [ circular: parseFinalVerdictCircular(ev) ], ev ] }

        // Non-zero OATK exits are propagated by the process itself so Nextflow's
        // retry/error strategy remains effective. This channel therefore contains
        // only structured successful outcomes: assembled or a valid no-contig result.
        ch_oatk_status = OATK.out.status

        ch_versions = ch_versions.mix(OATK.out.versions.first())
        ch_versions = ch_versions.mix(OATK_CHECK.out.versions.first())
        ch_summary_files = ch_summary_files.mix(OATK.out.log.map { _meta, log -> log })
        ch_summary_files = ch_summary_files.mix(OATK.out.fasta.map { _meta, fasta -> fasta })
        ch_summary_files = ch_summary_files.mix(OATK.out.gfa.map { _meta, gfa -> gfa })
        ch_summary_files = ch_summary_files.mix(ch_oatk_status.map { _meta, status -> status })
        // Oatk check evidence feeds the assembly summary (circularised override +
        // anomaly/length review reason), same as the GetOrganelle/MitoHiFi check evidence.
        ch_summary_files = ch_summary_files.mix(OATK_CHECK.out.evidence.map { _meta, ev -> ev })
        ch_multiqc_files = ch_multiqc_files.mix(OATK_CHECK.out.tool_params.collect { it[1] })
    } else {
        ASSEMBLY_NO_RESULT(ch_direct_oatk_reads.map { meta, _reads, reason -> [meta, reason] })
        ch_routed_failure_fasta = ASSEMBLY_NO_RESULT.out.fasta
        ch_routed_failure_log = ASSEMBLY_NO_RESULT.out.status
        ch_summary_files = ch_summary_files.mix(ASSEMBLY_NO_RESULT.out.status.map { _meta, status -> status })
        ch_versions = ch_versions.mix(ASSEMBLY_NO_RESULT.out.versions.first())
    }

    //
    // Fold the circularity verdict into meta.circular.
    //   The verdict comes from the check-circularity evidence (final_verdict_circular);
    //   failed assemblies have no evidence and carry circular=null (unknown). Every
    //   emitted channel that is later joined by the whole meta map downstream
    //   (assembly_fasta <-> assembly_log in the parent workflow; circularity_evidence
    //   in UPLOAD_RESULTS; reference_gb in MITOGENOME_ANNOTATION) is enriched from the
    //   SAME verdict so those joins still match.
    //
    ch_circ_verdict = MITOHIFI_CHECK_CIRCULARITY.out.evidence
        .map { meta, tsv -> [ meta, parseFinalVerdictCircular(tsv) ] }

    ch_assembly_fasta = MITOHIFI_MITOHIFI.out.fasta.mix(ch_routed_failure_fasta)
        .join(ch_circ_verdict, by: 0, remainder: true)
        .map { meta, fasta, circ -> [ meta + [ circular: circ ], fasta ] }

    ch_assembly_log = ch_assembly_log.mix(ch_routed_failure_log)
        .join(ch_circ_verdict, by: 0, remainder: true)
        .map { meta, log, circ -> [ meta + [ circular: circ ], log ] }

    ch_reference_gb = RELABEL_REFERENCE_GB.out.gb
        .join(ch_circ_verdict, by: 0, remainder: true)
        .map { meta, gb, circ -> [ meta + [ circular: circ ], gb ] }

    ch_circularity_evidence = MITOHIFI_CHECK_CIRCULARITY.out.evidence
        .map { meta, evidence -> [ meta + [ circular: parseFinalVerdictCircular(evidence) ], evidence ] }

    //
    // Collect MultiQC inputs and versions
    //   - Prefer feeding human‑readable logs and summary tables to MultiQC.
    //

    // MitoHiFi per-sample stats and logs (circularity-corrected stats)
    ch_multiqc_files = ch_multiqc_files.mix(ch_assembled_stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_CHECK_CIRCULARITY.out.evidence.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_CHECK_CIRCULARITY.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_AVERAGE_COVERAGE.out.coverage.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_MITOHIFI.out.logs.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_MITOHIFI.out.command_logs.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_FINDMITOREFERENCE.out.tool_params.collect { it[1] })
    ch_summary_files = ch_summary_files.mix(MITOHIFI_FINDMITOREFERENCE.out.status.map { _meta, status -> status })
    ch_multiqc_files = ch_multiqc_files.mix(CAT_FASTQ.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_MITOHIFI.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MITOHIFI_AVERAGE_COVERAGE.out.tool_params.collect { it[1] })
    ch_summary_files = ch_summary_files.mix(MITOHIFI_MITOHIFI.out.fasta.map { meta, fasta -> fasta })
    ch_summary_files = ch_summary_files.mix(MITOHIFI_MITOHIFI.out.gb.map { meta, gb -> gb })
    ch_summary_files = ch_summary_files.mix(MITOHIFI_MITOHIFI.out.command_logs.map { meta, log -> log })
    ch_summary_files = ch_summary_files.mix(MITOHIFI_MITOHIFI.out.logs.map { meta, log -> log })
    // The findMitoReference GenBank carries the reference species + accession the
    // assembly summary reports. It cannot be staged under its native name (the
    // accession-named `.gb` inside an identically-named `MitoReference` dir collides
    // across samples in MITOGENOME_ASSEMBLY_SUMMARY's flat collect()), so relabel it
    // to a per-sample `<assembly_prefix>.reference.gb` first.
    ch_summary_files = ch_summary_files.mix(RELABEL_REFERENCE_GB.out.gb.map { _meta, gb -> gb })
    // Pre-assembly reference-divergence flag (<prefix>.reference_divergence.txt):
    // the summary strips the suffix to the run prefix and folds a non-congeneric
    // verdict into the manual_review_reason.
    ch_summary_files = ch_summary_files.mix(REFERENCE_DIVERGENCE.out.flag.map { _meta, flag -> flag })
    // Re-selection audit trail: which candidates were considered and why one won,
    // plus the candidate-lookup outcome for samples that produced none.
    ch_summary_files = ch_summary_files.mix(REFERENCE_RANK.out.ranking.map { _meta, ranking -> ranking })
    ch_summary_files = ch_summary_files.mix(REFERENCE_CANDIDATES.out.status.map { _meta, status -> status })
    ch_summary_files = ch_summary_files.mix(ch_assembled_stats.map { meta, stats -> stats })
    // The circularity-check evidence is a per-run sidecar: the assembly summary
    // strips its .circularity_check.tsv suffix to the run prefix (so it joins the
    // existing run rather than spawning a phantom) and folds its length/repeat
    // anomaly into the manual_review_reason.
    ch_summary_files = ch_summary_files.mix(MITOHIFI_CHECK_CIRCULARITY.out.evidence.map { meta, evidence -> evidence })
    ch_summary_files = ch_summary_files.mix(MITOHIFI_AVERAGE_COVERAGE.out.coverage.map { meta, coverage -> coverage })

    // Versions for versions.yml collation (not MultiQC inputs)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    ch_versions = ch_versions.mix(MITOHIFI_MITOHIFI.out.versions.first())
    ch_versions = ch_versions.mix(MITOHIFI_AVERAGE_COVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(MITOHIFI_CHECK_CIRCULARITY.out.versions.first())
    ch_versions = ch_versions.mix(REFERENCE_DIVERGENCE.out.versions.first())
    ch_versions = ch_versions.mix(REFERENCE_CANDIDATES.out.versions.first())
    ch_versions = ch_versions.mix(REFERENCE_RANK.out.versions.first())

    // Per-sample bundle of everything this stage publishes into <prefix>/mtdna,
    // keyed by the (original) assembly prefix so the collapse mirror can restage the
    // whole folder into <prefix>_collapsed/mtdna for genuinely collapsed samples.
    // Excludes MITOHIFI_AVERAGE_COVERAGE.out.stats: it shares a basename with the
    // corrected CHECK_CIRCULARITY stats, which overwrites it in the real mtdna dir.
    ch_mtdna_files = MITOHIFI_MITOHIFI.out.fasta
        .mix( MITOHIFI_MITOHIFI.out.stats,
              MITOHIFI_MITOHIFI.out.gb,
              MITOHIFI_MITOHIFI.out.logs,
              MITOHIFI_MITOHIFI.out.command_logs,
              MITOHIFI_AVERAGE_COVERAGE.out.coverage,
              MITOHIFI_CHECK_CIRCULARITY.out.stats,
              MITOHIFI_CHECK_CIRCULARITY.out.evidence,
              REFERENCE_DIVERGENCE.out.flag )
        .map { meta, f -> [ meta.mt_assembly_prefix, f ] }
        .groupTuple()

    // Fold the oatk fallback's mtdna files (assembly + circularity check) into the same
    // bundle, keyed by the oatk assembly prefix, so a genuinely collapsed oatk concatemer
    // gets its full provenance mirrored into <prefix>_collapsed/mtdna like the other
    // assemblers. Empty (no oatk fallback ran, or the fallback is disabled).
    ch_oatk_mtdna_files = ch_oatk_fasta
        .mix(ch_oatk_circularity_evidence)
        .map { meta, f -> [ meta.mt_assembly_prefix, f ] }
        .groupTuple()
    ch_mtdna_files = ch_mtdna_files.mix(ch_oatk_mtdna_files)


    //
    // Emit outputs
    //

    emit:
    mtdna_files     = ch_mtdna_files               // channel: [ mt_assembly_prefix, [ mtdna files ] ]
    assembly_fasta  = ch_assembly_fasta            // channel: [ meta(+circular), assembly.fasta ]
    oatk_fasta      = ch_oatk_fasta                // channel: [ meta(+circular), oatk.mito.ctg.fasta ] (empty unless fallback enabled)
    oatk_log        = ch_oatk_log                  // channel: [ meta, oatk.log ]
    assembly_log    = ch_assembly_log              // channel: [ meta(+circular), contigs_stats.tsv ]
    reference_gb    = ch_reference_gb              // channel: [ meta(+circular), reference.gb ]
    // Full merged HiFi reads keyed by assembly prefix, for MITOGENOME_COVERAGE in the
    // parent workflow. Keyed rather than meta-joined because meta gains `circular`
    // downstream, so a whole-meta join would never match. Oatk is mixed in from its own
    // channel: ch_oatk_input rewrites mt_assembly_prefix, so its reads are not reachable
    // through ch_reads_by_prefix.
    depth_reads     = ch_reads_by_prefix.map { m, r -> [ m.mt_assembly_prefix, r ] }
                        .mix(ch_oatk_reads.map { m, r -> [ m.mt_assembly_prefix, r ] })
    // Fold the oatk fallback's circularity evidence into the same channel so the parent
    // workflow's collapse-concatemer join, QC gate and assembly summary treat oatk exactly
    // like MitoHiFi (matched join item -> flows to annotation incrementally). Empty unless
    // the fallback ran.
    circularity_evidence = ch_circularity_evidence.mix(ch_oatk_circularity_evidence) // channel: [ meta(+circular), *_check.tsv ]
    summary_files   = ch_summary_files
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
