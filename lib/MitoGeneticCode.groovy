/*
 * The NCBI mitochondrial genetic code (translation table) for a taxonomic class.
 *
 * Mitochondrial codes differ by lineage and the wrong one is not a cosmetic error: it
 * mistranslates every CDS in MITOS2 annotation and QC protein translation, and it is
 * written as the mgcode into the table2asn submission, so it reaches ENA/GenBank.
 *
 * This lookup lives in lib/ rather than inside a subworkflow because two entrypoints
 * need it and they resolve the class from different places: the main pipeline takes it
 * from the enriched samplesheet, while qc_only_from_annotations.nf has no samplesheet
 * and queries it from SQL. Both must reach the same code for the same class, otherwise
 * re-QCing an assembly through the standalone entrypoint would silently rewrite its
 * mgcode.
 *
 * Only the class -> code lookup is shared. What to do when a class is NOT mapped is
 * deliberately left to the caller, because the two callers differ: the main pipeline
 * aborts on an invertebrate with no confirmed code (it is about to produce the
 * annotation, so a wrong table is unrecoverable), while the QC-only entrypoint falls
 * back to --translation_table with a warning (it is re-QCing annotations that already
 * exist, and refusing to run helps nobody).
 */
class MitoGeneticCode {

    // Cnidaria: the Coelenterate code.
    static final Set<String> CODE_4_CLASSES = [
        'anthozoa', 'hydrozoa', 'scyphozoa', 'cubozoa', 'staurozoa',
        'myxozoa', 'polypodiozoa',
    ] as Set

    // Echinoderms and flatworms: the Echinoderm/Flatworm code.
    static final Set<String> CODE_9_CLASSES = [
        'asteroidea', 'ophiuroidea', 'echinoidea', 'holothuroidea', 'crinoidea',
        'rhabditophora', 'trematoda', 'cestoda', 'monogenea', 'turbellaria',
    ] as Set

    private static String norm(taxClass) {
        return (taxClass == null) ? '' : taxClass.toString().trim().toLowerCase()
    }

    /*
     * The confirmed code for a class, or null when the class is unmapped or unresolved.
     *
     * null means "no confirmed code", never "use the default": vertebrates and
     * not-yet-curated invertebrate lineages are indistinguishable here, and only the
     * caller knows whether defaulting is safe in its context.
     */
    static Integer forClass(taxClass) {
        def c = norm(taxClass)
        if (!c || c == 'unknown') return null
        if (c in CODE_4_CLASSES) return 4
        if (c in CODE_9_CLASSES) return 9
        return null
    }
}
