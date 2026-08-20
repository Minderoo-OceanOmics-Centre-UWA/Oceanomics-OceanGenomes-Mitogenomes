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
// Aliased so conf/modules.config can address this invocation on its own. The same
// module also runs in MITOGENOME_ANNOTATION_LCA, on the PUBLISHED assembly (post
// collapse/curation) -- that one is the curator-facing verdict and keeps the publish
// rule. This one grades the RAW MitoHiFi output, early enough to route on, and is not
// published: it is a routing signal like countGenbankCds, and publishing it would put a
// second, differently-scoped reference_relevance.txt in the sample's annotation dir.
include { REFERENCE_RELEVANCE as REFERENCE_RELEVANCE_ROUTING } from '../../../../modules/local/reference_relevance'
include { OATK                             } from '../../../../modules/local/oatk'
include { OATK_CHECK                       } from '../../../../modules/local/oatk/check_circularity'
include { ASSEMBLY_NO_RESULT               } from '../../../../modules/local/assembly_no_result'

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

// Did REFERENCE_RANK actually choose a reference? It emits EMPTY chosen_reference
// files when it declined -- no candidate recruited any reads, no reads were
// subsampled, or a lone candidate went unmapped -- because a substitution made on no
// read evidence cannot be known to improve on the reference findMitoReference already
// resolved. OG56 is the case in point: assembled against a Fowleria vaiulae record
// picked out of a five-way all-zero tie, at 76.2% identity to its own assembly.
// Truthiness alone is not enough here; the files exist either way.
def hasChosenReference(fasta, gb) {
    return fasta && gb && fasta.size() > 0 && gb.size() > 0
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
//
// ROUTING SIGNAL ONLY. The single caller (the gene-incomplete oatk branch below) uses
// this as a filter predicate and then drops the value: it is never published, never
// written to a file, and never reported in the assembly summary, MultiQC or the SQL
// upload.
//
// It is NOT the pipeline's protein-coding-gene count. These are CDS features written
// by MitoHiFi's own reference-guided annotation, whereas every assembly is re-annotated
// downstream by EMMA / MITOS2; the reported num_cds comes from that re-annotation via
// *.annotation_stats.csv (see bin/mitogenome_assembly_summary.py). The two numbers can
// legitimately disagree and neither overwrites the other -- they only share the
// params.mitogenome_summary_expected_pcg_count threshold. Reading MitoHiFi's GenBank is
// the right source *for this decision*: it measures what reference-guided read
// recruitment lost, which is exactly the collapse being detected.
def countGenbankCds(gb) {
    try {
        if (!gb || !gb.exists() || gb.size() == 0) return null
        return gb.text.readLines().count { it.startsWith('     CDS ') }
    } catch (ignored) {
        return null
    }
}

// Read one named column from the single data row of a per-sample TSV. Returns null
// when the file, the column or the cell is missing, which every caller treats as "no
// evidence" rather than as a defect -- an unreadable artefact must never divert an
// assembly on its own.
def firstRowColumn(tsv, name) {
    try {
        if (!tsv || !tsv.exists() || tsv.size() == 0) return null
        def rows = tsv.text.readLines().findAll { it?.trim() }
        if (rows.size() < 2) return null
        def idx = rows[0].split('\t', -1).findIndexOf { it.trim() == name }
        if (idx < 0) return null
        def cells = rows[1].split('\t', -1)
        return idx < cells.size() ? cells[idx].trim() : null
    } catch (ignored) {
        return null
    }
}

// Assembly length as a multiple of the reference length, from MITOHIFI_CHECK_CIRCULARITY.
// A reference too distant to recruit cleanly can inflate an assembly as readily as it
// can truncate one: OG56 came back circular with all 13 CDS at 1.281x reference length,
// carrying a ~170 bp control-region unit in ~78.6 tandem copies. The CDS count sees
// nothing wrong with that molecule, so length is the second, independent measure of the
// same failure. Null when unparseable (see firstRowColumn).
//
// ROUTING SIGNAL ONLY, exactly as countGenbankCds above: read as a predicate and
// discarded. The value a curator reads comes from the published circularity_check.tsv
// this parses, which already reports length_ratio, excess_bp and the repeat geometry.
def parseLengthRatio(tsv) {
    try {
        def raw = firstRowColumn(tsv, 'length_ratio')
        return raw ? raw as Double : null
    } catch (ignored) {
        return null
    }
}

// PASS / DIVERGENT / MISMATCH / UNKNOWN from a REFERENCE_RELEVANCE flag file, which
// grades the resolved reference against the assembly it produced. Same one-line
// tab-separated shape as REFERENCE_DIVERGENCE, so it parses the same way.
def parseRelevanceVerdict(txt) {
    try {
        return txt.text.readLines().find { it?.trim() }?.split('\t', -1)?.first()?.trim()?.toUpperCase() ?: 'UNKNOWN'
    } catch (ignored) {
        return 'UNKNOWN'
    }
}

// Decide whether a MitoHiFi assembly is defective enough to warrant the reference-free
// Oatk fallback, and say which defect caught it. Returns null to leave the assembly
// alone. Named and file-scoped so tests/assembly_routing can exercise the real decision
// rather than a copy of it (see that harness for why a copied predicate is worse than
// no test at all).
//
// Judges the ASSEMBLY, not the reference that produced it. Two independent defects, from
// artefacts MitoHiFi and the circularity check already write:
//
//   gene_incomplete  Recruitment against a too-distant reference drops the most divergent
//                    gene blocks and returns a clean-looking but truncated molecule.
//                    OG2102 (Rouleina attrita, non-congeneric Alepocephalus reference)
//                    came back as a plausible 15.6 kb contig carrying 7 of 13 PCGs, so it
//                    never reached the empty-FASTA fallback.
//
//   length_inflated  The same distant reference can inflate instead of truncate. OG56
//                    (Vincentia punctata, Fowleria reference at 76.2% identity) came back
//                    circular with all 13 CDS at 1.281x reference length, carrying a
//                    ~170 bp control-region unit in ~78.6 tandem copies. The CDS count
//                    sees nothing wrong with it, which is why length is checked separately.
//
// length_inflated is QUALIFIED by the relevance verdict and gene_incomplete is NOT, and the
// asymmetry is load-bearing. Over-length is ambiguous on its own: OG848, OG852 and OG853 all
// run 1.12-1.33x against congeneric references that PASS, which is real length heteroplasmy
// and not something a different assembler will "fix". A PASS is therefore evidence the length
// is the organism's, not the reference's. Missing genes admit no such reading -- OG2102 is
// PASS at 7 of 13 PCGs -- so qualifying gene_incomplete would silently un-route the very
// sample the branch was built for.
//
// Every condition demands POSITIVE evidence: a missing or unreadable artefact yields null
// and the assembly is left alone. That is also why the qualifier tests for DIVERGENT/MISMATCH
// rather than `!= PASS` -- UNKNOWN means the check reached no verdict, which is not evidence
// the reference caused anything.
//
// Order is precedence, not preference: OG56 satisfies more than one condition and the caller
// stamps every routed sample with the same oatk prefix, so a sample that produced two reasons
// would launch two OATK tasks writing identical output names. First match wins makes that
// unrepresentable.
def oatkFallbackReason(gb, evidence, relevance, expected_pcg_count, length_ratio_threshold) {
    def cds     = countGenbankCds(gb)
    def ratio   = parseLengthRatio(evidence)
    def verdict = relevance ? parseRelevanceVerdict(relevance) : 'UNKNOWN'

    if (cds != null && cds < expected_pcg_count) {
        return 'gene_incomplete_mitohifi'
    }
    if (ratio != null && ratio >= length_ratio_threshold && verdict in ['DIVERGENT', 'MISMATCH']) {
        return 'length_inflated_mitohifi'
    }
    if (verdict == 'MISMATCH') {
        // The reference neither covers nor matches the assembly, i.e. it is simply the
        // wrong record -- the reference-free case by definition. Untested on real data:
        // the only MISMATCH observed so far (OG810) is a GetOrganelle sample, which never
        // reaches this branch.
        return 'wrong_reference_mitohifi'
    }
    return null
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

    // Stand-in for assemblies that never reached a circularity check (failed / no-contig).
    // Keeps circularity_evidence one-row-per-assembly; see the emit block for why totality
    // is a contract rather than a convenience.
    def no_circularity_evidence = file("${projectDir}/assets/empty_circularity_check.tsv", checkIfExists: true)

    //
    // Canonical read channel. The assembly prefix is embedded once, up front, and
    // both the reference lookup and every read route reuse that same meta map. The
    // prefix is part of the meta, so it is part of every downstream join key:
    // deriving it a second time (from a meta that a process may have handed back)
    // risks the two copies drifting apart and silently emptying a join. One
    // derivation, one key.
    //
    // The version is parsed from the pinned container tag (params.mitohifi_container)
    // rather than `mitohifi.py --version`, because the 3.2.3 release ships a stale
    // self-reported version (3.2.1). The tag is the single source of truth: bump it in
    // nextflow.config and the version in every assembly name follows automatically.
    //

    // Two fields with different lifetimes -- see the GetOrganelle subworkflow for the full
    // rationale. mt_assembly_prefix is IDENTITY (re-stamped to the FASTA basename wherever
    // curation renames the assembly, e.g. _collapsed / _concat); mt_assembly_run_prefix is
    // LINEAGE, written once here and never reassigned, for the joins that reunite artefacts
    // which are per sample+assembler and identical across variants.
    //
    // Note this stamp is EARLIER in its subworkflow than GetOrganelle's: the reference routing
    // below already keys on the prefix, so CAT_FASTQ and MITOHIFI_FINDMITOREFERENCE sit
    // downstream of it. Do not move it later to reclaim their cache -- that would be a second
    // key-lifetime change layered on this one.
    ch_reads_prefixed = fastp_reads
        .map { meta, reads ->
            def mt_assembly_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${mitohifi_version_stripped}mitohifi"
            [ meta + [
                mt_assembly_prefix:     mt_assembly_prefix,
                mt_assembly_run_prefix: mt_assembly_prefix
            ], reads ]
        }

    // map just the meta for the species query
    ch_species_reference = ch_reads_prefixed.map { meta, _reads -> meta }
    // .view()

    //
    // MODULE: Find a closely related species for reference
    //

    MITOHIFI_FINDMITOREFERENCE (
        ch_species_reference
    )

    //
    // Route reads on whether they actually need concatenating.
    //
    // Only groups that carry more than one FASTQ (more than one pair, when paired)
    // are worth a CAT_FASTQ task; a group holding a single file is already the merged
    // FASTQ MitoHiFi wants. `passthrough` is an unconditional fallback rather than the
    // negation of the first condition, so the two branches are mutually exclusive and
    // jointly exhaustive by construction: no sample can fall through the branch and
    // disappear, which is how the singleton route was lost once before.
    //

    ch_reads_routed = ch_reads_prefixed.branch { meta, reads ->
        def readList = reads instanceof List ? reads : [ reads ]
        needs_concat: meta.single_end ? readList.size() > 1 : readList.size() > 2
        passthrough: true
    }

    //
    // MODULE: Concatenate fastq reads where there are multiple fastq files
    //

    CAT_FASTQ (
        ch_reads_routed.needs_concat
    )

    // Recombine the two routes with `mix`. `mix` forwards each item the moment it
    // arrives, so singleton samples reach the reference join, MitoHiFi and Oatk while
    // the multi-file groups are still concatenating. Anything that has to see a whole
    // channel first (collect / groupTuple / a remainder join) would instead hold every
    // sample back until the slowest concatenation finished.
    final_reads = ch_reads_routed.passthrough.mix(CAT_FASTQ.out.reads)

    //
    // Combine fastp files with the mito reference output
    //
    // Keyed on mt_assembly_prefix, NOT on the whole meta map. On a resumed run
    // findMitoReference's meta comes back from the cache database, and a cache-restored
    // meta is not guaranteed to `equals` the live one that never went through a process:
    // Nextflow's task hash ignores an empty collection, so an optional samplesheet column
    // the sheet does not carry (nf-schema materialises it as `[]`) changes meta shape
    // without changing the hash. The task still resumes, hands back the older meta, and a
    // whole-meta join silently matches nothing -- which is exactly how every singleton
    // HiFi sample stopped reaching MitoHiFi. The prefix is one stable string per assembly,
    // and the reads-side meta (live, canonical, already carrying the prefix) is what
    // travels downstream. .toString() because the prefix is a GString and a GString never
    // equals a String.
    //
    // Reads keep the meta they were routed with, so they stay joinable even when
    // findMitoReference emits a no-reference placeholder.
    //

    ch_reference_outcomes = MITOHIFI_FINDMITOREFERENCE.out.reference
        .join(MITOHIFI_FINDMITOREFERENCE.out.status, by: 0)
        .map { meta, ref_fasta, ref_gb, status ->
            [ meta.mt_assembly_run_prefix.toString(), ref_fasta, ref_gb, status ]
        }
    ch_reference_joined = final_reads
        .map { meta, reads -> [ meta.mt_assembly_run_prefix.toString(), meta, reads ] }
        .join(ch_reference_outcomes, by: 0)
        .map { _prefix, meta, reads, ref_fasta, ref_gb, status ->
            [ meta, reads, ref_fasta, ref_gb, status ]
        }
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
            hasChosenReference(chosen_fasta, chosen_gb)
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
    // final_reads because ch_oatk_input below OVERWRITES mt_assembly_prefix
    // with the oatk prefix, so a prefix-keyed join against final_reads would
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
        // Oatk assembler tag for the run prefix, parsed from the container version
        // (1.0 -> v10oatk) so the summary files it as an Oatk run distinct from the
        // failed MitoHiFi attempt. Overwrite mt_assembly_prefix so every downstream
        // output (Oatk contig, annotation, summary) is named and grouped under it.
        def oatk_version_stripped = (params.oatk_container ?: 'oatk:1.0')
            .tokenize(':').last().tokenize('--').first().replaceAll('\\.', '')

        // Reads keyed by the assembly prefix for the two joins below, for the same reason
        // the reference join is keyed that way: MITOHIFI_MITOHIFI's meta is cache-restored
        // on a resume while final_reads is computed live, and the two are not guaranteed
        // to `equals` even when they describe the same sample.
        ch_reads_keyed = final_reads.map { meta, reads -> [ meta.mt_assembly_run_prefix.toString(), reads ] }

        ch_failed_oatk_reads = ch_mitohifi_fasta_branched.failed
            .map { meta, _empty_fasta -> [ meta.mt_assembly_run_prefix.toString(), meta ] }
            .join(ch_reads_keyed, by: 0)
            .map { _prefix, meta, reads -> [meta, reads, 'empty_mitohifi'] }

        // An empty FASTA is not the only way a divergent reference ruins an assembly, and a
        // defective assembly is what the fallback is really keyed on -- so judge the ASSEMBLY,
        // not the reference that produced it. oatkFallbackReason (defined at file scope, with
        // the full rationale) reads the defect from artefacts MitoHiFi and the circularity
        // check already write, and names which one caught the sample.
        //
        // Every signal here is a filter predicate and is discarded; none is published. The
        // whole branch sits inside `if (params.enable_oatk_fallback)`, false by default, so on
        // a default run none of it is computed. Failure modes stay cheap in both directions: a
        // concatemer duplicates CDS features so the count lands above the threshold and never
        // trips the `<` filter, and a spurious route costs one extra Oatk run, because the
        // MitoHiFi assembly is not withdrawn (see below).
        def expected_pcg_count = (params.mitogenome_summary_expected_pcg_count ?: 13) as int
        def length_ratio_threshold = (params.mitogenome_oatk_length_ratio_threshold ?: 1.15) as Double

        // Grade the resolved reference against the assembly it produced, here rather than
        // waiting for MITOGENOME_ANNOTATION_LCA's copy: the verdict is needed to route, and
        // routing happens before annotation. No cycle is introduced -- the module needs only
        // a FASTA and the reference GenBank, both of which exist at this point.
        //
        // Invoked inside the enable_oatk_fallback gate so a default run pays nothing for a
        // signal it would not act on.
        //
        // Keyed on the run prefix, not the whole meta: RELABEL_REFERENCE_GB and
        // MITOHIFI_MITOHIFI are different processes whose restored metas need not `equals`
        // on a -resume, the failure documented for the oatk reference join below. The
        // assembly's meta is the one carried forward, so out.flag joins cleanly against the
        // other MitoHiFi outputs. Inner join is total here: every assembled sample came
        // through ch_reference_resolved, which is also what RELABEL_REFERENCE_GB consumes.
        ch_relabelled_gb_keyed = RELABEL_REFERENCE_GB.out.gb
            .map { m, gb -> [ m.mt_assembly_run_prefix.toString(), gb ] }

        REFERENCE_RELEVANCE_ROUTING (
            ch_mitohifi_fasta_branched.assembled
                .map { meta, fasta -> [ meta.mt_assembly_run_prefix.toString(), meta, fasta ] }
                .join(ch_relabelled_gb_keyed, by: 0)
                .map { _prefix, meta, fasta, gb -> [ meta, fasta, gb ] }
        )

        // One channel, one reason per sample. ch_oatk_input below stamps every entry with the
        // SAME oatk_prefix, so a sample arriving on two channels would launch two OATK tasks
        // writing identical output names; resolving the reason before the mix makes that
        // unrepresentable, which is why this is not three mixed channels.
        //
        // The relevance join carries remainder:true while the GenBank and evidence joins do
        // not: REFERENCE_RELEVANCE_ROUTING runs errorStrategy 'ignore', and an inner join
        // would let one failed blastn silently withdraw gene_incomplete routing from a sample
        // like OG2102, which does not depend on the relevance verdict at all.
        ch_defective_oatk_reads = ch_mitohifi_fasta_branched.assembled
            .join(MITOHIFI_MITOHIFI.out.gb, by: 0)
            .join(MITOHIFI_CHECK_CIRCULARITY.out.evidence, by: 0)
            .join(REFERENCE_RELEVANCE_ROUTING.out.flag, by: 0, remainder: true)
            .filter { it[1] != null }   // drop any right-only remainder
            .map { items ->
                [ items[0], oatkFallbackReason(
                    items[2],                                  // MitoHiFi GenBank
                    items[3],                                  // circularity evidence
                    items.size() > 4 ? items[4] : null,        // relevance flag, may be absent
                    expected_pcg_count,
                    length_ratio_threshold) ]
            }
            .filter { _meta, reason -> reason != null }
            .map { meta, reason -> [ meta.mt_assembly_run_prefix.toString(), meta, reason ] }
            .join(ch_reads_keyed, by: 0)
            .map { _prefix, meta, reason, reads -> [meta, reads, reason] }

        // The MitoHiFi assembly is deliberately NOT withdrawn when this fires: it stays
        // published alongside the Oatk attempt (distinct v10oatk prefix), so the summary
        // shows both and curation picks. Discarding it would lose information.
        //
        // fallback_reason records WHY oatk was invoked (no_reference /
        // reference_lookup_error / cross_order_reference / empty_mitohifi /
        // gene_incomplete_mitohifi / length_inflated_mitohifi / wrong_reference_mitohifi).
        // It is carried in meta for provenance only: nothing in the repo reads it, so the
        // routing decision currently leaves no trace in any output file. Surfacing it (and
        // the CDS count / length ratio behind it) as summary evidence is a deliberate TODO,
        // not an oversight to be fixed in passing -- and it matters more now that the
        // post-assembly gate has three ways to fire and a first-match precedence between
        // them, so which one caught a sample is no longer inferable from the outputs.
        // Oatk is a distinct assembly RUN, not a curated variant of the MitoHiFi one: it has
        // its own assembler prefix, its own reads entry and its own DB row. So it re-stamps
        // the lineage key too, rather than inheriting MitoHiFi's -- otherwise its artefacts
        // would join against the MitoHiFi run they were meant to replace.
        ch_oatk_input = ch_direct_oatk_reads
            .mix(ch_failed_oatk_reads, ch_defective_oatk_reads)
            .map { meta, reads, reason ->
                def oatk_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${oatk_version_stripped}oatk"
                [ meta + [
                    mt_assembly_prefix:     oatk_prefix,
                    mt_assembly_run_prefix: oatk_prefix,
                    circular: null,
                    assembler_fallback: 'oatk',
                    fallback_reason: reason
                ], reads ]
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
        // circularity read. Keyed on mt_assembly_prefix, not the whole meta:
        // OATK.out.fasta and OATK.out.gfa both carry cache-restored meta on a
        // -resume, and a whole-meta join drops every sample whose restored meta no
        // longer `equals` the live one -- the same failure the reference join above
        // documents, and exactly why OG109/OG2089 never reached OATK_CHECK.
        ch_oatk_fa_gfa = ch_oatk_assembled
            .map { m, fasta -> [ m.mt_assembly_prefix.toString(), m, fasta ] }
            .join(OATK.out.gfa.map { m, gfa -> [ m.mt_assembly_prefix.toString(), gfa ] }, by: 0)
            .map { _prefix, m, fasta, gfa -> [ m, fasta, gfa ] }   // [meta, fasta, gfa]

        // Attach the per-sample reference GenBank. RELABEL_REFERENCE_GB.out.gb carries the
        // MitoHiFi-prefixed meta (not the oatk one), so key on [id, sequencing_type, date].
        // The true no-reference oatk path (ch_direct_oatk_reads: ch_reference_branched.missing
        // + the cross-order route) has no relabelled reference, so pair every such sample with
        // an explicit NO_REFERENCE placeholder (check_getorganelle.py treats a zero-length
        // reference as absent). RELABEL runs on ch_reference_resolved = keep + reselect, which
        // is disjoint from ch_direct_oatk_reads, so the two contributors never double-key a
        // sample; together they make ch_oatk_ref_keyed TOTAL over the assembled oatk contigs
        // and the join below becomes a plain per-sample join (was a whole-run remainder wait).
        def no_reference_gb_oatk = file("${projectDir}/assets/NO_REFERENCE.gb", checkIfExists: true)
        ch_oatk_ref_keyed = RELABEL_REFERENCE_GB.out.gb
            .map { m, gb -> [ [m.id, m.sequencing_type, m.date], gb ] }
            .mix( ch_direct_oatk_reads
                    .map { meta, _reads, _reason -> [ [meta.id, meta.sequencing_type, meta.date], no_reference_gb_oatk ] } )
        ch_oatk_check_in = ch_oatk_fa_gfa
            .map { m, fasta, gfa -> [ [m.id, m.sequencing_type, m.date], m, fasta, gfa ] }
            .join(ch_oatk_ref_keyed, by: 0)
            .map { _key, m, fasta, gfa, ref -> [ m, fasta, gfa, ref ] }   // [meta, fasta, gfa, ref]

        OATK_CHECK ( ch_oatk_check_in )

        // Fold the corrected circular verdict into meta on every oatk channel that is later
        // joined by the whole meta map (fasta <-> evidence in the parent collapse join; the
        // QC-gate evidence). Keyed on mt_assembly_prefix, not the whole meta: OATK.out.* and
        // OATK_CHECK.out.evidence carry meta from different cache entries on a -resume, and a
        // whole-meta join drops the verdict for every sample whose restored meta no longer
        // `equals` the live one (same failure the reference join above documents). remainder:true
        // keeps the no-contig oatk (no evidence -> circular:null, empty FASTA filtered out of
        // annotation downstream).
        // One verdict row per emitted oatk FASTA. OATK_CHECK runs on assembled contigs only
        // (ch_oatk_fa_gfa), so the no-contig empties carry a null verdict, emitted here from the
        // same local emptiness filter. That makes ch_oatk_circ_verdict TOTAL over OATK.out.fasta
        // (and over OATK.out.log below), so both consumers use a plain per-sample join instead
        // of a remainder join that waits on the whole run to close.
        ch_oatk_circ_verdict = OATK_CHECK.out.evidence
            .map { m, tsv -> [ m.mt_assembly_prefix.toString(), parseFinalVerdictCircular(tsv) ] }
            .mix( OATK.out.fasta
                    .filter { _m, fasta -> fasta.size() == 0 }
                    .map { m, _fasta -> [ m.mt_assembly_prefix.toString(), null ] } )

        ch_oatk_fasta = OATK.out.fasta
            .map { m, fasta -> [ m.mt_assembly_prefix.toString(), m, fasta ] }
            .join(ch_oatk_circ_verdict, by: 0)
            .map { _prefix, m, fasta, circ -> [ m + [ circular: circ ], fasta ] }

        // Push input for oatk: the getorg_check.tsv evidence when OATK_CHECK ran (so the push
        // script parses final_verdict_circular into a real stats value), else the raw .oatk.log
        // so a no-contig / no-evidence oatk still emits a tuple and hits the push script's
        // fail-loud path rather than vanishing. NOT OATK.out.log alone, which is the raw syncasm
        // log that matches no parser and silently wrote stats=NULL. The raw log stays in
        // ch_summary_files (below) for the summary / MultiQC.
        // Evidence keyed per oatk sample, made TOTAL over OATK.out.log by pairing the no-contig
        // empties (no OATK_CHECK evidence) with a null row so `ev ?: log` falls back to the raw
        // oatk.log exactly as the old remainder path did -- but with a plain join that releases
        // each sample immediately.
        ch_oatk_evidence_keyed = OATK_CHECK.out.evidence
            .map { m, ev -> [ m.mt_assembly_prefix.toString(), ev ] }
            .mix( OATK.out.fasta
                    .filter { _m, fasta -> fasta.size() == 0 }
                    .map { m, _fasta -> [ m.mt_assembly_prefix.toString(), null ] } )
        ch_oatk_push_input = OATK.out.log
            .map { m, log -> [ m.mt_assembly_prefix.toString(), m, log ] }
            .join(ch_oatk_evidence_keyed, by: 0)
            .map { prefix, m, log, ev -> [ prefix, m, ev ?: log ] }
        ch_oatk_log = ch_oatk_push_input
            .join(ch_oatk_circ_verdict, by: 0)
            .map { _prefix, m, push_file, circ -> [ m + [ circular: circ ], push_file ] }

        // One evidence row per emitted oatk FASTA, never fewer. OATK_CHECK runs on the
        // assembled contigs only (ch_oatk_assembled), so the no-contig empties are paired
        // with the empty-check placeholder here instead of being left absent. Downstream
        // joins the evidence per sample rather than waiting for the channel to close, so a
        // sample with no evidence row would be dropped at the QC gate rather than delayed.
        // The complement is taken from the same local filter that built ch_oatk_assembled,
        // so no join (and no channel-close dependency) is needed to identify it.
        ch_oatk_circularity_evidence = OATK_CHECK.out.evidence
            .map { m, ev -> [ m + [ circular: parseFinalVerdictCircular(ev) ], ev ] }
            .mix( OATK.out.fasta
                    .filter { _meta, fasta -> fasta.size() == 0 }
                    .map { m, _fasta -> [ m + [ circular: null ], no_circularity_evidence ] } )

        // Non-zero OATK exits are propagated by the process itself so Nextflow's
        // retry/error strategy remains effective. This channel therefore contains
        // only structured successful outcomes: assembled or a valid no-contig result.
        ch_oatk_status = OATK.out.status

        ch_versions = ch_versions.mix(OATK.out.versions.first())
        ch_versions = ch_versions.mix(OATK_CHECK.out.versions.first())
        ch_versions = ch_versions.mix(REFERENCE_RELEVANCE_ROUTING.out.versions.first())
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
    // One verdict row per assembly this arm emits, never fewer. MITOHIFI_CHECK_CIRCULARITY
    // runs on the `assembled` branch only, so the `failed` branch and the routed failures
    // (ASSEMBLY_NO_RESULT) carry no evidence -- pair them here with a null verdict (unknown,
    // not "not circular"). Both complements come from branches/channels already local at this
    // point, so establishing totality costs no join and introduces no channel-close dependency.
    // With ch_circ_verdict total over every left key below, those joins become plain per-sample
    // joins that release each finished sample immediately instead of waiting for the whole run
    // to close (the remainder-join barrier this replaces).
    ch_circ_verdict = MITOHIFI_CHECK_CIRCULARITY.out.evidence
        .map { meta, tsv -> [ meta, parseFinalVerdictCircular(tsv) ] }
        .mix( ch_mitohifi_fasta_branched.failed.map { meta, _fasta -> [ meta, null ] } )
        .mix( ch_routed_failure_fasta.map { meta, _fasta -> [ meta, null ] } )

    ch_assembly_fasta = MITOHIFI_MITOHIFI.out.fasta.mix(ch_routed_failure_fasta)
        .join(ch_circ_verdict, by: 0)
        .map { meta, fasta, circ -> [ meta + [ circular: circ ], fasta ] }

    ch_assembly_log = ch_assembly_log.mix(ch_routed_failure_log)
        .join(ch_circ_verdict, by: 0)
        .map { meta, log, circ -> [ meta + [ circular: circ ], log ] }

    ch_reference_gb = RELABEL_REFERENCE_GB.out.gb
        .join(ch_circ_verdict, by: 0)
        .map { meta, gb, circ -> [ meta + [ circular: circ ], gb ] }

    // One evidence row per emitted assembly FASTA, never fewer. MITOHIFI_CHECK_CIRCULARITY
    // runs on the `assembled` branch only, and ch_assembly_fasta additionally carries the
    // `failed` branch plus the routed failures, so those two are paired with the empty-check
    // placeholder here (circular stays null -- unknown, not "not circular"). Both complements
    // come from branches/channels that are already local at this point, so establishing
    // totality costs no join and introduces no channel-close dependency.
    //
    // Totality matters downstream: the parent workflow and the QC gate now attach evidence
    // with a plain per-key join instead of collecting the whole channel, which is what lets a
    // finished sample cross the gate while other assemblies are still running. A missing
    // evidence row would silently drop its sample there rather than merely delay it.
    ch_circularity_evidence = MITOHIFI_CHECK_CIRCULARITY.out.evidence
        .map { meta, evidence -> [ meta + [ circular: parseFinalVerdictCircular(evidence) ], evidence ] }
        .mix( ch_mitohifi_fasta_branched.failed
                .map { meta, _fasta -> [ meta + [ circular: null ], no_circularity_evidence ] } )
        .mix( ch_routed_failure_fasta
                .map { meta, _fasta -> [ meta + [ circular: null ], no_circularity_evidence ] } )

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
    // Sized with groupKey so an assembled sample's bundle releases as soon as its own files
    // arrive instead of waiting for the whole run to close. A fully assembled (therefore
    // reference-guided) sample contributes exactly nine files -- fasta, stats, gb, hifiasm log,
    // command log, average coverage, corrected stats, circularity evidence and the divergence
    // flag -- one per mixed channel, which is the maximum any sample reaches. groupKey emits
    // each such bundle THE MOMENT it reaches nine, before the channel closes, so every bundle
    // the collapse mirror can actually consume streams -- the collapse mirror only ever
    // restages ASSEMBLED samples (a failed empty assembly has no molecule to collapse).
    //
    // This bundle's count is genuinely non-uniform (unlike the GetOrganelle/oatk bundles and
    // the ENA record, which are exactly fixed): a FAILED assembly contributes fewer files (no
    // coverage/check, no gb). Those failed bundles are never consumed downstream, but rather
    // than drop them (a plain groupKey would, and a rare assembled sample missing one optional
    // output would be dropped with it -- silent provenance loss), remainder:true flushes every
    // short group at close, exactly as the bare groupTuple did. So this remainder is NOT the
    // 2e anti-pattern: no consumed bundle is delayed by it -- assembled bundles already streamed
    // at nine -- it only backstops the unused failed ones and guards against loss. Nine is the
    // ceiling, so no sample overshoots and no bundle is ever split.
    ch_mtdna_files = MITOHIFI_MITOHIFI.out.fasta
        .mix( MITOHIFI_MITOHIFI.out.stats,
              MITOHIFI_MITOHIFI.out.gb,
              MITOHIFI_MITOHIFI.out.logs,
              MITOHIFI_MITOHIFI.out.command_logs,
              MITOHIFI_AVERAGE_COVERAGE.out.coverage,
              MITOHIFI_CHECK_CIRCULARITY.out.stats,
              MITOHIFI_CHECK_CIRCULARITY.out.evidence,
              REFERENCE_DIVERGENCE.out.flag )
        .map { meta, f -> [ groupKey(meta.mt_assembly_run_prefix, 9), f ] }
        .groupTuple(remainder: true)
        .map { key, files -> [ key.getGroupTarget(), files ] }

    // Fold the oatk fallback's mtdna files (assembly + circularity check) into the same
    // bundle, keyed by the oatk assembly prefix, so a genuinely collapsed oatk concatemer
    // gets its full provenance mirrored into <prefix>_collapsed/mtdna like the other
    // assemblers. Empty (no oatk fallback ran, or the fallback is disabled).
    // Exactly two files per oatk sample (the contig + its circularity evidence, both now
    // total over every oatk sample), so groupKey(prefix, 2) releases each bundle as soon as
    // both arrive rather than at whole-run close. The count is uniform, so a plain groupTuple
    // releases every sample and no close-time remainder is needed. Unwrap the GroupKey so the
    // collapse mirror's key join still matches.
    ch_oatk_mtdna_files = ch_oatk_fasta
        .mix(ch_oatk_circularity_evidence)
        .map { meta, f -> [ groupKey(meta.mt_assembly_run_prefix, 2), f ] }
        .groupTuple()
        .map { key, files -> [ key.getGroupTarget(), files ] }
    ch_mtdna_files = ch_mtdna_files.mix(ch_oatk_mtdna_files)


    //
    // Emit outputs
    //

    emit:
    mtdna_files     = ch_mtdna_files               // channel: [ mt_assembly_run_prefix, [ mtdna files ] ]
    assembly_fasta  = ch_assembly_fasta            // channel: [ meta(+circular), assembly.fasta ]
    oatk_fasta      = ch_oatk_fasta                // channel: [ meta(+circular), oatk.mito.ctg.fasta ] (empty unless fallback enabled)
    oatk_log        = ch_oatk_log                  // channel: [ meta(+circular), getorg_check.tsv (assembled) | oatk.log (no-contig) ] — push input
    assembly_log    = ch_assembly_log              // channel: [ meta(+circular), contigs_stats.tsv ]
    reference_gb    = ch_reference_gb              // channel: [ meta(+circular), reference.gb ]
    // The complete HiFi read set per sample (concatenated where the group held more than one
    // FASTQ, otherwise the single original file), keyed by the LINEAGE prefix, for
    // MITOGENOME_COVERAGE in the parent workflow. Keyed rather than meta-joined because meta
    // gains `circular` and a re-stamped identity downstream, so a whole-meta join would never
    // match; lineage rather than identity because these are reads, not an assembly, and there
    // is one set per sample+assembler whatever curation later renames the molecule to. Oatk is
    // mixed in from its own channel: it is a separate assembly RUN with its own lineage key, so
    // its reads are not reachable through final_reads.
    depth_reads     = final_reads.map { m, r -> [ m.mt_assembly_run_prefix, r ] }
                        .mix(ch_oatk_reads.map { m, r -> [ m.mt_assembly_run_prefix, r ] })
    // Fold the oatk fallback's circularity evidence into the same channel so the parent
    // workflow's collapse-concatemer join, QC gate and assembly summary treat oatk exactly
    // like MitoHiFi (matched join item -> flows to annotation incrementally). Empty unless
    // the fallback ran.
    //
    // CONTRACT: exactly one row per assembly FASTA this subworkflow emits (assembly_fasta +
    // oatk_fasta), using assets/empty_circularity_check.tsv where no check ran. The parent
    // workflow and the QC gate rely on that totality to attach evidence with a plain per-key
    // join; a channel that is merely "evidence where it exists" would force them back to a
    // remainder join or a collect, which is what previously pinned every finished sample to
    // the slowest assembly in the run.
    circularity_evidence = ch_circularity_evidence.mix(ch_oatk_circularity_evidence) // channel: [ meta(+circular), *_check.tsv ]
    summary_files   = ch_summary_files
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
