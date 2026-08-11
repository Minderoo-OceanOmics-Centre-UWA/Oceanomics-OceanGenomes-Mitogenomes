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
// canonical row carried the final variant's prefix. Real runs key the canonical row on the
// SAMPLE-level prefix -- verified against run mitogenomes-missing-audit-5, where OG868
// produced exactly two rows: canonical "OG868.hic.250624.getorg1770" carrying the reseed
// FASTA and a real depth, plus provenance "OG868.hic.250624.getorg1770reseed" with the
// placeholder. The duplicate also could not run at all, so nobody noticed.)
include { selectProvenanceVariants } from '../../subworkflows/local/mitogenome_assembly/getorganelle/main.nf'

workflow UPLOAD_SELECTION {
    take:
    variant_inputs // [ variant_prefix, meta, fasta, log ]
    checked_circ   // [ sample_prefix, checked_variant_prefix, verdict ]

    main:
    results = selectProvenanceVariants(variant_inputs, checked_circ)

    emit:
    results
}

workflow {
    if (params.scenario == 'upload_selection') {
        // meta.mt_assembly_prefix is the SAMPLE-level prefix on every variant row; the
        // reseed / _rgj suffix lives only in the variant prefix (the FASTA basename).
        variants = Channel.of(
            ['OG1.getorg1770',          [id: 'OG1', mt_assembly_prefix: 'OG1.getorg1770'], 'raw-first.fa',  'first.log'],
            ['OG1.getorg1770reseed',    [id: 'OG1', mt_assembly_prefix: 'OG1.getorg1770'], 'raw-reseed.fa', 'reseed.log'],
            ['OG1.getorg1770reseed_rgj',[id: 'OG1', mt_assembly_prefix: 'OG1.getorg1770'], 'raw-rgj.fa',    'rgj.log']
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
