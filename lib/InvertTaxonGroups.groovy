/**
 * Taxonomic class groupings used to route invertebrate samples (meta.invertebrates
 * == true) through phylum-appropriate annotation logic, instead of the
 * originally coral-only (Cnidaria) machinery.
 *
 * Two groups matter, and they are NOT the same set:
 *   - REDUCED_TRNA: Cnidaria + Porifera. Both groups legitimately have
 *     reduced/atypical mitochondrial tRNA complements (corals lack tRNA-Phe;
 *     sponge tRNA counts range from 2 to 27 across lineages, with documented
 *     tRNA-Phe loss in some clades). They share the Coelenterate/Mold genetic
 *     code (4) and need cox1-anchored re-origin instead of the usual
 *     tRNA-Phe-based anchor.
 *   - CORAL_FIX_ELIGIBLE: Cnidaria only. The ANNOTATION_QC_GATE /
 *     CORAL_ANNOTATION_FIX pathway repairs a specific failure mode (dropped
 *     16S, nad5 split across the giant nad5-717 group I intron) that is a
 *     Hexacorallia-specific trait, absent in Octocorallia and not shared by
 *     Porifera (whose own group I introns, when present, sit in cox1, in only
 *     a few sponge families, via lineage-specific horizontal transfer -- a
 *     different gene and a different mechanism). Porifera therefore stays out
 *     of this group even though it is in REDUCED_TRNA above.
 */
class InvertTaxonGroups {

    static final Set<String> CNIDARIA_CLASSES = [
        'anthozoa', 'hydrozoa', 'scyphozoa', 'cubozoa', 'staurozoa', 'myxozoa', 'polypodiozoa',
    ] as Set

    static final Set<String> PORIFERA_CLASSES = [
        'demospongiae', 'calcarea', 'hexactinellida', 'homoscleromorpha', 'porifera',
    ] as Set

    static final Set<String> ECHINODERM_FLATWORM_CLASSES = [
        'asteroidea', 'ophiuroidea', 'echinoidea', 'holothuroidea', 'crinoidea',
        'rhabditophora', 'trematoda', 'cestoda', 'monogenea', 'turbellaria',
    ] as Set

    // Echinoderm-only subset of the above, used for cox1 rotation-panel selection
    // (flatworms are not represented in the curated panels yet and fall back to
    // the generic anchor -- ROTATE_ORIGIN no-ops safely on a weak hit).
    static final Set<String> ECHINODERMATA_CLASSES = [
        'asteroidea', 'ophiuroidea', 'echinoidea', 'holothuroidea', 'crinoidea', 'echinodermata',
    ] as Set

    static final Set<String> MOLLUSCA_CLASSES = [
        'bivalvia', 'gastropoda', 'cephalopoda', 'polyplacophora', 'scaphopoda',
        'monoplacophora', 'caudofoveata', 'solenogastres', 'mollusca',
    ] as Set

    static final Set<String> ARTHROPODA_CLASSES = [
        'malacostraca', 'hexanauplia', 'branchiopoda', 'ostracoda',
        'cephalocarida', 'remipedia', 'maxillopoda', 'pycnogonida', 'arthropoda',
    ] as Set

    static final Set<String> REDUCED_TRNA_CLASSES = CNIDARIA_CLASSES + PORIFERA_CLASSES

    static final Set<String> CORAL_FIX_ELIGIBLE_CLASSES = CNIDARIA_CLASSES

    private static String norm(taxClass) {
        (taxClass ?: '').toString().trim().toLowerCase()
    }

    static boolean isReducedTrna(taxClass) {
        norm(taxClass) in REDUCED_TRNA_CLASSES
    }

    static boolean isCoralFixEligible(taxClass) {
        norm(taxClass) in CORAL_FIX_ELIGIBLE_CLASSES
    }

    // Which curated cox1 anchor panel (assets/cox1_<group>.faa) a sample's
    // ROTATE_ORIGIN re-origin step should tblastn against. 'reduced_trna'
    // (Cnidaria/Porifera) is also the fallback for any invertebrate class not in
    // one of the curated bilaterian panels -- rotate_to_cox1.py passes the
    // assembly through unrotated on a weak/no hit, so a mismatched panel is a
    // safe no-op, never a corruption.
    static String cox1PanelGroup(taxClass) {
        def c = norm(taxClass)
        if (c in REDUCED_TRNA_CLASSES) return 'reduced_trna'
        if (c in MOLLUSCA_CLASSES) return 'mollusca'
        if (c in ARTHROPODA_CLASSES) return 'arthropoda'
        if (c in ECHINODERMATA_CLASSES) return 'echinodermata'
        return 'reduced_trna'
    }
}
