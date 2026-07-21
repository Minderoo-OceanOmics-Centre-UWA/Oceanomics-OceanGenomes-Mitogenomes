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

workflow UPLOAD_SELECTION {
    take:
    raw_variants
    canonical_results

    main:
    raw_candidates = raw_variants.map { meta, fasta, log ->
        [ meta.mt_assembly_prefix, [ priority: 0, meta: meta, fasta: fasta, log: log ] ]
    }
    canonical_candidates = canonical_results.map { meta, fasta, log ->
        [ meta.mt_assembly_prefix, [ priority: 1, meta: meta, fasta: fasta, log: log ] ]
    }
    selected = raw_candidates
        .mix(canonical_candidates)
        .groupTuple(by: 0)
        .map { _prefix, candidates ->
            def result = candidates.max { it.priority }
            [ result.meta, result.fasta, result.log ]
        }

    emit:
    results = selected
}

workflow {
    if (params.scenario == 'upload_selection') {
        raw = Channel.of(
            [[id: 'OG1', mt_assembly_prefix: 'OG1.first'], 'raw-first.fa', 'first.log'],
            [[id: 'OG1', mt_assembly_prefix: 'OG1.reseed'], 'raw-reseed.fa', 'reseed.log'],
            [[id: 'OG1', mt_assembly_prefix: 'OG1.rgj'], 'raw-rgj.fa', 'rgj.log']
        )
        canonical = Channel.of(
            [[id: 'OG1', mt_assembly_prefix: 'OG1.rgj', circular: true], 'curated-rgj.fa', 'rgj.log'],
            [[id: 'OG2', mt_assembly_prefix: 'OG2.oatk'], 'oatk.fa', 'oatk.log']
        )
        UPLOAD_SELECTION(raw, canonical)
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
