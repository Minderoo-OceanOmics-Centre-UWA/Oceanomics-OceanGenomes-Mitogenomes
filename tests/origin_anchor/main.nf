nextflow.enable.dsl = 2

// Exercises the taxon -> published-origin anchor lookup that the annotation
// subworkflow passes to MITOS2 and CORAL_ANNOTATION_FIX as --origin-gene.
//
// The Anthozoa rows are the point of this harness. Before the anchor table every
// invertebrate was re-origined to tRNA-Met unconditionally; the table replaces that
// with a measurement, and the measurement only works if it is keyed on ORDER, because
// Anthozoa as a class has no majority convention and its four orders each have a
// different one. A class-only lookup would rotate every submitted stony coral off
// tRNA-Met and change its ENA sequence checksum.
workflow ORIGIN_ANCHOR {
    take:
    taxa   // channel of [ sample_id, tax_order, tax_class ]

    main:
    InvertTaxonGroups.loadOriginAnchors(
        file("${projectDir}/assets/taxonomy/mito_origin_anchors.json", checkIfExists: true))

    resolved = taxa.map { sample_id, tax_order, tax_class ->
        tuple(sample_id, tax_order, tax_class,
              InvertTaxonGroups.originAnchor(tax_order, tax_class))
    }

    emit:
    resolved
}

workflow {
    ORIGIN_ANCHOR(Channel.of(
        ['S1', 'Scleractinia', 'Anthozoa'],
        ['S2', 'Malacalcyonacea', 'Anthozoa'],
        ['S3', '', 'Demospongiae'],
        ['S4', 'Valvatida', 'Asteroidea']))
    ORIGIN_ANCHOR.out.resolved.view { id, o, c, a -> "ANCHOR\t${id}\t${o}\t${c}\t${a}" }
}
