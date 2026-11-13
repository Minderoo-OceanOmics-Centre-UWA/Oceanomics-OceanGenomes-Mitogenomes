/*
 * The NCBI mitochondrial genetic code (translation table) for a taxonomic class.
 *
 * Mitochondrial codes differ by lineage and the wrong one is not a cosmetic error: it
 * mistranslates every CDS in MITOS2 annotation and QC protein translation, and it is
 * written as the mgcode into the table2asn submission, so it reaches ENA/GenBank.
 *
 * This is the front door two entrypoints share, because they resolve the class from
 * different places: the main pipeline takes it from the enriched samplesheet, while
 * qc_only_from_annotations.nf has no samplesheet and queries it from SQL. Both must
 * reach the same code for the same class, otherwise re-QCing an assembly through the
 * standalone entrypoint would silently rewrite its mgcode.
 *
 * The map itself lives in assets/taxonomy/mito_genetic_codes.json and is owned by
 * InvertTaxonGroups, so bin/create_samplesheet.py can resolve the same codes when it
 * writes the samplesheet's genetic_code column. This class adds nothing to it but the
 * loaded-ness check below.
 *
 * Only the lookup is shared. What to do when a class is NOT mapped is deliberately left
 * to the caller, because the two callers differ: the main pipeline aborts on an
 * invertebrate with no confirmed code (it is about to produce the annotation, so a wrong
 * table is unrecoverable), while the QC-only entrypoint falls back to --translation_table
 * with a warning (it is re-QCing annotations that already exist, and refusing to run
 * helps nobody).
 */
class MitoGeneticCode {

    /*
     * Parse the class -> code asset. Idempotent, so a caller can invoke it without a
     * guard, and every entrypoint that resolves a code must call it first.
     */
    static void load(codesFile) {
        InvertTaxonGroups.loadGeneticCodes(codesFile)
    }

    /*
     * The confirmed code for a class, or null when the class is unmapped or unresolved.
     *
     * null means "no confirmed code", never "use the default": vertebrates, the
     * documented ambiguous classes, and not-yet-curated lineages are indistinguishable
     * here, and only the caller knows whether defaulting is safe in its context.
     *
     * An unloaded map is an error rather than a null. InvertTaxonGroups.geneticCode()
     * returns null when the asset was never parsed, which is indistinguishable from a
     * genuinely unmapped class -- an entrypoint that forgot to load would silently give
     * every sample the caller's fallback, which on the QC-only path is exactly the
     * blanket vertebrate code that path was fixed to stop using.
     */
    static Integer forClass(taxClass) {
        if (!InvertTaxonGroups.geneticCodesLoaded()) {
            throw new IllegalStateException(
                "MitoGeneticCode.forClass() called before the class -> code map was " +
                "parsed -- call MitoGeneticCode.load(" +
                "file(\"\${projectDir}/assets/taxonomy/mito_genetic_codes.json\")) first")
        }
        def c = (taxClass == null) ? '' : taxClass.toString().trim()
        // 'unknown' is the literal bin/create_samplesheet.py writes for an unresolved
        // taxon, and it must not be looked up as if it were a class name.
        if (!c || c.equalsIgnoreCase('unknown')) return null
        return InvertTaxonGroups.geneticCode(c)
    }
}
