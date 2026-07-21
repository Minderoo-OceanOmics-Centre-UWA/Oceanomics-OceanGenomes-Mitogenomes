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
include { OATK                             } from '../../../../modules/local/oatk'
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
    // MODULE: Relabel the reference GenBank to a per-sample name so it can be fed
    // (collision-free) to the assembly summary, which reads the reference species
    // and accession from it.
    //

    RELABEL_REFERENCE_GB (
        combined_with_mt_assembly_prefix.map { meta, _fasta, _ref_fasta, ref_gb -> [meta, ref_gb] }
    )

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

    REFERENCE_DIVERGENCE (
        RELABEL_REFERENCE_GB.out.gb
    )

    ch_reference_route = combined_with_mt_assembly_prefix
        .join(REFERENCE_DIVERGENCE.out.flag, by: 0)
        .branch { _meta, _reads, _ref_fasta, _ref_gb, flag ->
            oatk_direct: parseReferenceTier(flag) == 'CROSS_ORDER'
            mitohifi: true
        }


    //
    // MODULE: Run assembly using MitoHifi from reads
    //

    MITOHIFI_MITOHIFI (
        ch_reference_route.mitohifi.map { meta, reads, ref_fasta, ref_gb, _flag -> [meta, reads, ref_fasta, ref_gb] },
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
    //   Runs only on samples MitoHiFi failed to assemble; emits a FASTA that the
    //   parent workflow feeds into annotation. Off by default (needs an Oatk
    //   container + an OatkDB mito profile for the sample clade).
    //
    ch_oatk_fasta = Channel.empty()
    ch_oatk_log = Channel.empty()
    ch_routed_failure_fasta = Channel.empty()
    ch_routed_failure_log = Channel.empty()

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

        ch_oatk_input = ch_direct_oatk_reads.mix(ch_failed_oatk_reads)
            .map { meta, reads, reason ->
                def oatk_prefix = "${meta.id}.${meta.sequencing_type}.${meta.date}.v${oatk_version_stripped}oatk"
                [ meta + [ mt_assembly_prefix: oatk_prefix, circular: null, assembler_fallback: 'oatk', fallback_reason: reason ], reads ]
            }

        // Stage the .fam AND its .h3* nhmmer indexes together (nhmmscan needs them
        // side by side); params.oatk_mito_db points at the .fam.
        ch_oatk_db = Channel.fromPath("${params.oatk_mito_db}*", checkIfExists: true).collect()

        OATK (
            ch_oatk_input,
            ch_oatk_db
        )

        ch_oatk_fasta = OATK.out.fasta
        ch_oatk_log = OATK.out.log
        // Non-zero OATK exits are propagated by the process itself so Nextflow's
        // retry/error strategy remains effective. This channel therefore contains
        // only structured successful outcomes: assembled or a valid no-contig result.
        ch_oatk_status = OATK.out.status

        ch_versions = ch_versions.mix(OATK.out.versions.first())
        ch_summary_files = ch_summary_files.mix(OATK.out.log.map { _meta, log -> log })
        ch_summary_files = ch_summary_files.mix(OATK.out.fasta.map { _meta, fasta -> fasta })
        ch_summary_files = ch_summary_files.mix(OATK.out.gfa.map { _meta, gfa -> gfa })
        ch_summary_files = ch_summary_files.mix(ch_oatk_status.map { _meta, status -> status })
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


    //
    // Emit outputs
    //

    emit:
    assembly_fasta  = ch_assembly_fasta            // channel: [ meta(+circular), assembly.fasta ]
    oatk_fasta      = ch_oatk_fasta                // channel: [ meta(+circular), oatk.mito.ctg.fasta ] (empty unless fallback enabled)
    oatk_log        = ch_oatk_log                  // channel: [ meta, oatk.log ]
    assembly_log    = ch_assembly_log              // channel: [ meta(+circular), contigs_stats.tsv ]
    coverage_stats  = MITOHIFI_AVERAGE_COVERAGE.out.coverage
    reference_gb    = ch_reference_gb              // channel: [ meta(+circular), reference.gb ]
    circularity_evidence = ch_circularity_evidence // channel: [ meta(+circular), circularity_check.tsv ]
    summary_files   = ch_summary_files
    multiqc_files   = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions        = ch_versions              // channel: [ path(versions.yml) ]
}
