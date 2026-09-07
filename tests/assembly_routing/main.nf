nextflow.enable.dsl = 2

include { MITOHIFI_FINDMITOREFERENCE } from '../../modules/nf-core/mitohifi/findmitoreference/main'

workflow REFERENCE_POLICY {
    take:
    samples

    main:
    MITOHIFI_FINDMITOREFERENCE(samples)

    emit:
    reference = MITOHIFI_FINDMITOREFERENCE.out.reference
    status = MITOHIFI_FINDMITOREFERENCE.out.status
}

// Import the real implementation rather than copying it, so this cannot silently pass
// against a stale duplicate of the logic it is meant to protect. (It previously WAS such a
// duplicate: a copy of the old groupTuple selection, fed hand-written inputs in which the
// canonical row carried the final variant's prefix. The duplicate also could not run at all,
// so nobody noticed.)
//
// The canonical row now carries the WINNING assembly's own name, so this selection excludes
// the winner and keeps every superseded attempt. It used to exclude the first pass, which was
// correct only while the canonical row was misnamed -- filed under the sample-level prefix
// even when the molecule in it was a reseed. That misnaming is what dropped curated assemblies
// at the QC gate and filed their depth against a superseded row.
include { selectProvenanceVariants } from '../../subworkflows/local/mitogenome_assembly/getorganelle/main.nf'

// Same reasoning as above: import the real Oatk-fallback predicate, not a copy. It decides
// whether a MitoHiFi assembly is defective enough to re-assemble reference-free, from files
// on disk, so the tests hand it real (tiny) artefacts rather than pre-parsed values -- that
// way the parsers it depends on (countGenbankCds, parseLengthRatio, parseRelevanceVerdict)
// are exercised too, and a change to the TSV shapes cannot pass silently.
include { oatkFallbackReason } from '../../subworkflows/local/mitogenome_assembly/mitohifi/main.nf'

// Same reasoning again: import the real partial-failure predicate and the real
// was_circular reader, not copies. These two decide whether a MitoHiFi assembly that
// crashed midway is carried through the subworkflow at all, and what topology it is
// recorded with. OG2133 was lost because the first of them did not exist.
include { mitohifiPartialFailure; mitohifiStatsCircular } from '../../subworkflows/local/mitogenome_assembly/mitohifi/main.nf'

workflow MITOHIFI_PARTIAL_FAILURE {
    take:
    cases   // [ label, command_log ]

    main:
    results = cases.map { label, log -> [ label, mitohifiPartialFailure(log) ] }

    emit:
    results
}

workflow MITOHIFI_STATS_CIRCULAR {
    take:
    cases   // [ label, contigs_stats ]

    main:
    // Groovy null does not survive the channel round-trip as a distinguishable value in
    // the assertions, so map it to the literal 'null' string here. The distinction that
    // matters is null-vs-false, and both are visible this way.
    results = cases.map { label, stats ->
        def v = mitohifiStatsCircular(stats)
        [ label, v == null ? 'null' : v.toString() ]
    }

    emit:
    results
}

workflow OATK_FALLBACK_ROUTING {
    take:
    cases   // [ label, gb, evidence, relevance ]

    main:
    results = cases.map { label, gb, evidence, relevance ->
        [ label, oatkFallbackReason(gb, evidence, relevance, 13, 1.15) ]
    }

    emit:
    results
}

workflow UPLOAD_SELECTION {
    take:
    variant_inputs // [ variant_prefix(identity), meta(+mt_assembly_run_prefix), fasta, log ]
    checked_circ   // [ mt_assembly_run_prefix, checked_identity, verdict ]

    main:
    results = selectProvenanceVariants(variant_inputs, checked_circ)

    emit:
    results
}

workflow {
    if (params.scenario == 'upload_selection') {
        // meta.mt_assembly_run_prefix is the LINEAGE key, identical on every variant row; each
        // variant's own identity is the key (the FASTA basename).
        variants = Channel.of(
            ['OG1.getorg1770',          [id: 'OG1', mt_assembly_run_prefix: 'OG1.getorg1770'], 'raw-first.fa',  'first.log'],
            ['OG1.getorg1770reseed',    [id: 'OG1', mt_assembly_run_prefix: 'OG1.getorg1770'], 'raw-reseed.fa', 'reseed.log'],
            ['OG1.getorg1770reseed_rgj',[id: 'OG1', mt_assembly_run_prefix: 'OG1.getorg1770'], 'raw-rgj.fa',    'rgj.log']
        )
        checked = Channel.of(
            ['OG1.getorg1770', 'OG1.getorg1770reseed_rgj', true]
        )
        UPLOAD_SELECTION(variants, checked)
        UPLOAD_SELECTION.out.results.view { meta, fasta, log ->
            "RESULT\t${meta.mt_assembly_prefix}\t${fasta}\t${meta.circular}\t${log}"
        }
    } else {
        def species = params.scenario ?: 'success'
        REFERENCE_POLICY(Channel.value([
            id: species,
            species_id: species,
            nominal_species_id: species,
            mt_assembly_prefix: species
        ]))
        REFERENCE_POLICY.out.status.view { meta, status -> "STATUS\t${meta.id}\t${status.text.readLines()[1]}" }
    }
}
