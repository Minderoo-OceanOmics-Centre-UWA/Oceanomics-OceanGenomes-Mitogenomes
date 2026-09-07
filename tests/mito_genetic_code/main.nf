nextflow.enable.dsl = 2

// Exercises the class -> mitochondrial genetic code lookup that both the main
// pipeline (via mitoGeneticCode() in prepare_samplesheet) and the QC-only
// entrypoint resolve through. The two used to disagree by construction: the
// QC-only path had no lookup at all and gave every sample --translation_table.
workflow MITO_GENETIC_CODE {
    take:
    classes   // channel of [ sample_id, tax_class ]

    main:
    MitoGeneticCode.load(
        file("${projectDir}/assets/taxonomy/mito_genetic_codes.json", checkIfExists: true))

    resolved = classes.map { sample_id, tax_class ->
        tuple(sample_id, tax_class, MitoGeneticCode.forClass(tax_class))
    }

    emit:
    resolved
}

workflow {
    MITO_GENETIC_CODE(Channel.of(['OG1', 'Anthozoa'], ['OG2', 'Asteroidea'], ['OG3', 'Actinopteri'],
                                 ['OG4', 'Bivalvia'], ['OG5', 'Ascidiacea'], ['OG6', 'Pterobranchia']))
    MITO_GENETIC_CODE.out.resolved.view { id, c, code -> "CODE\t${id}\t${c}\t${code}" }
}
