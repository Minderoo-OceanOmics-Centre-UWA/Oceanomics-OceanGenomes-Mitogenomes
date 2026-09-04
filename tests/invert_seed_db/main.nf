// Harness for InvertTaxonGroups.seedDbGroup(): resolve a class to its curated
// GetOrganelle seed database and report both the group and whether the database
// files that GETORGANELLE_RESEED would stage actually exist on disk.
//
// InvertTaxonGroups lives in lib/, so it is on the classpath here exactly as it
// is in the real subworkflow -- this cannot pass against a stale copy.

nextflow.enable.dsl = 2

workflow INVERT_SEED_DB {
    take:
    tax_classes   // channel of class names

    main:
    resolved = tax_classes.map { taxClass ->
        def group = InvertTaxonGroups.seedDbGroup(taxClass)
        def seed  = group ? file("${projectDir}/assets/refdb/${group}/${group}_mito_refdb.fasta")       : null
        def genes = group ? file("${projectDir}/assets/refdb/${group}/${group}_mito_refdb.label.fasta") : null
        [ taxClass.toString().trim(), group ?: 'NONE',
          (seed && seed.exists()) ? 'seed:yes' : 'seed:no',
          (genes && genes.exists()) ? 'genes:yes' : 'genes:no' ].join(' ')
    }

    emit:
    resolved
}
