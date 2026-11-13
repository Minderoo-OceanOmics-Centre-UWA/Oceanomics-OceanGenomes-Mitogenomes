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
 *
 * Mitochondrial genetic codes are deliberately NOT one of these sets. They live
 * in assets/taxonomy/mito_genetic_codes.json, loaded here by loadGeneticCodes(), because
 * bin/create_samplesheet.py needs the same mapping and a second copy in Python
 * would drift. They also do not line up with the sets above -- the code-4 group
 * includes Ctenophora, which is neither reduced-tRNA nor coral-fix eligible, and
 * the flatworm classes split across codes 9, 14 and 21 -- so deriving one from
 * the other would corrupt panel selection and coral-fix routing.
 */
import groovy.json.JsonSlurper
class InvertTaxonGroups {

    // Each set also carries its own phylum name, because a sample whose taxonomy
    // only resolved to phylum rank arrives with the phylum in meta.class -- the
    // same convention assets/taxonomy/mito_genetic_codes.json uses.
    static final Set<String> CNIDARIA_CLASSES = [
        'anthozoa', 'hydrozoa', 'scyphozoa', 'cubozoa', 'staurozoa', 'myxozoa', 'polypodiozoa',
        'cnidaria',
    ] as Set

    static final Set<String> PORIFERA_CLASSES = [
        'demospongiae', 'calcarea', 'hexactinellida', 'homoscleromorpha', 'porifera',
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

    // Thecostraca (barnacles), Copepoda, Ichthyostraca and Mystacocarida are the
    // classes the old Maxillopoda has been split into; NCBI returns them, and
    // without them a barnacle got no seed database and the coral cox1 panel.
    static final Set<String> ARTHROPODA_CLASSES = [
        'malacostraca', 'hexanauplia', 'branchiopoda', 'ostracoda',
        'cephalocarida', 'remipedia', 'maxillopoda', 'pycnogonida', 'arthropoda',
        'thecostraca', 'copepoda', 'ichthyostraca', 'mystacocarida',
    ] as Set

    static final Set<String> CTENOPHORA_CLASSES = [
        'tentaculata', 'nuda', 'ctenophora',
    ] as Set

    static final Set<String> TUNICATA_CLASSES = [
        'ascidiacea', 'thaliacea', 'appendicularia', 'tunicata',
    ] as Set

    static final Set<String> ANNELIDA_CLASSES = [
        'polychaeta', 'clitellata', 'annelida', 'echiura', 'sipuncula',
    ] as Set

    static final Set<String> REDUCED_TRNA_CLASSES = CNIDARIA_CLASSES + PORIFERA_CLASSES

    static final Set<String> CORAL_FIX_ELIGIBLE_CLASSES = CNIDARIA_CLASSES

    private static String norm(taxClass) {
        (taxClass ?: '').toString().trim().toLowerCase()
    }

    // ---- Mitochondrial genetic codes (assets/taxonomy/mito_genetic_codes.json) ----

    private static Map<String, Integer> geneticCodes = null
    private static Map<String, String> ambiguousCodes = null

    /**
     * Parse assets/taxonomy/mito_genetic_codes.json once per session. Idempotent, so the
     * per-sample closure in prepare_samplesheet can call it without a guard.
     * A class listed twice with different codes, or listed both as a code and as
     * ambiguous, is a contradiction in the asset and stops the run here rather
     * than resolving to whichever entry happened to be parsed last.
     */
    static synchronized void loadGeneticCodes(codesFile) {
        if (geneticCodes != null) return
        def handle = codesFile instanceof File ? codesFile : new File(codesFile.toString())
        def payload = new JsonSlurper().parse(handle)

        def resolved = [:]
        payload.codes.each { entry ->
            def code = entry.code as int
            entry.classes.each { taxClass ->
                def key = norm(taxClass)
                if (resolved.containsKey(key) && resolved[key] != code) {
                    throw new IllegalStateException(
                        "${handle}: class '${key}' is mapped to both genetic code " +
                        "${resolved[key]} and ${code}")
                }
                resolved[key] = code
            }
        }

        def ambiguous = [:]
        (payload.ambiguous ?: [:]).each { taxClass, reason ->
            ambiguous[norm(taxClass)] = reason.toString()
        }

        def clash = resolved.keySet().intersect(ambiguous.keySet())
        if (clash) {
            throw new IllegalStateException(
                "${handle}: class(es) ${clash.sort().join(', ')} are listed both " +
                "with a genetic code and as ambiguous")
        }

        geneticCodes = resolved
        ambiguousCodes = ambiguous
    }

    /**
     * Whether loadGeneticCodes() has run. geneticCode() returns null both for an
     * unmapped class and for a map that was never parsed, so a caller that needs to
     * tell those apart (MitoGeneticCode.forClass) asks here first.
     */
    static boolean geneticCodesLoaded() {
        geneticCodes != null
    }

    /** Resolved mitochondrial genetic code for a class, or null if unmapped. */
    static Integer geneticCode(taxClass) {
        geneticCodes?.get(norm(taxClass))
    }

    /**
     * Why a class is deliberately left unmapped, or null if it is not one of the
     * documented ambiguous cases. Used to explain the abort rather than just
     * reporting the class as unknown.
     */
    static String ambiguousReason(taxClass) {
        ambiguousCodes?.get(norm(taxClass))
    }

    static boolean isReducedTrna(taxClass) {
        norm(taxClass) in REDUCED_TRNA_CLASSES
    }

    static boolean isCoralFixEligible(taxClass) {
        norm(taxClass) in CORAL_FIX_ELIGIBLE_CLASSES
    }

    // Which curated cox1 anchor panel (assets/panels/cox1/<group>.faa) a sample's
    // ROTATE_ORIGIN re-origin step should tblastn against. 'reduced_trna'
    // (Cnidaria/Porifera) is also the fallback for any invertebrate class not in
    // one of the curated bilaterian panels -- rotate_to_cox1.py passes the
    // assembly through unrotated on a weak/no hit, so a mismatched panel is a
    // safe no-op, never a corruption. Note this catch-all is the opposite of
    // seedDbGroup()'s null below, and deliberately so: the cost of a wrong panel
    // is nothing, the cost of a wrong seed is a whole wasted reseed.
    static String cox1PanelGroup(taxClass) {
        def c = norm(taxClass)
        if (c in REDUCED_TRNA_CLASSES) return 'reduced_trna'
        if (c in MOLLUSCA_CLASSES) return 'mollusca'
        if (c in ARTHROPODA_CLASSES) return 'arthropoda'
        if (c in ECHINODERMATA_CLASSES) return 'echinodermata'
        return 'reduced_trna'
    }

    // ---- Published origin anchors (assets/taxonomy/mito_origin_anchors.json) ----
    //
    // Which gene the PUBLISHED mitogenome is re-origined to, passed to MITOS2 and
    // CORAL_ANNOTATION_FIX as --origin-gene. Distinct from cox1PanelGroup() above:
    // that picks the panel for the PRE-annotation rotation, whose only job is to move
    // the linearisation point off an intron-split gene so MITOS can annotate cleanly,
    // and cox1 is a fine universal choice for that. This is the POST-annotation
    // rotation that decides where the deposited sequence starts, and it is measured
    // per taxon rather than guessed.
    //
    // Generated, not curated -- the one difference from mito_genetic_codes.json.
    // bin/build_origin_anchor_table.py tallies which gene sits at position 1 across
    // every record in assets/refdb/ and writes the table wholesale.
    //
    // Keyed on ORDER first, then class, then group, because the class level is not
    // fine enough for Anthozoa: as a class it has no majority (rrnL 33.9%, cox1 28.5%,
    // trnM 17.6%), which hides four orders that each have a clear and DIFFERENT
    // convention -- Scleractinia trnM 63.8%, Malacalcyonacea rrnL 69.3%, Zoantharia
    // cox1 58.6%, Scleralcyonacea cox1 51.6%. Resolving stony corals from the class
    // aggregate would rotate every already-submitted one off tRNA-Met and change its
    // ENA sequence checksum for no reason.

    private static Map<String, String> originAnchorOrders = null
    private static Map<String, String> originAnchorClasses = null
    private static Map<String, String> originAnchorGroups = null
    private static String originAnchorDefault = null

    /**
     * Parse assets/taxonomy/mito_origin_anchors.json once per session. Idempotent, so a
     * per-sample closure can call it without a guard, exactly like loadGeneticCodes().
     * A duplicated key, or an anchor outside the MITOS feature vocabulary, stops the run
     * here rather than silently resolving to whichever entry parsed last.
     */
    static synchronized void loadOriginAnchors(anchorFile) {
        if (originAnchorOrders != null) return
        def handle = anchorFile instanceof File ? anchorFile : new File(anchorFile.toString())
        def payload = new JsonSlurper().parse(handle)

        def readLevel = { section, label ->
            def out = [:]
            (payload[section] ?: [:]).each { key, entry ->
                def k = norm(key)
                if (out.containsKey(k)) {
                    throw new IllegalStateException("${handle}: ${label} '${k}' listed twice")
                }
                def anchor = entry?.anchor?.toString()
                if (!isValidAnchor(anchor)) {
                    throw new IllegalStateException(
                        "${handle}: ${label} '${k}' anchor '${anchor}' is not a MITOS feature key")
                }
                out[k] = anchor
            }
            out
        }

        def orders  = readLevel('orders',  'order')
        def classes = readLevel('classes', 'class')
        def groups  = readLevel('groups',  'group')

        def fallback = payload.policy?.default_anchor?.toString()
        if (!isValidAnchor(fallback)) {
            throw new IllegalStateException(
                "${handle}: policy.default_anchor '${fallback}' is not a MITOS feature key")
        }

        // A table generated without --taxdump-dir has no order level, which would
        // resolve Scleractinia from the Anthozoa class aggregate and silently rotate
        // stony corals off trnM. That must not reach a run.
        def levels = (payload.policy?.levels ?: []).collect { it.toString() }
        if (!levels || levels[0] != 'order') {
            throw new IllegalStateException(
                "${handle}: policy.levels is ${levels} -- the ORDER level is missing, so " +
                "Scleractinia would resolve from the Anthozoa class aggregate (no majority) " +
                "and stony corals would be rotated off tRNA-Met. Regenerate the table with " +
                "bin/build_origin_anchor_table.py --taxdump-dir <taxdump>.")
        }

        originAnchorOrders = orders
        originAnchorClasses = classes
        originAnchorGroups = groups
        originAnchorDefault = fallback
    }

    /** Whether loadOriginAnchors() has run. */
    static boolean originAnchorsLoaded() {
        originAnchorOrders != null
    }

    // MITOS/EMMA feature keys an anchor is allowed to be. TL and TS are deliberately
    // absent: the reference databases record tRNA-Leu/tRNA-Ser without saying which
    // copy, so those tallies exist but are not addressable in a MITOS annotation.
    private static final Set<String> ANCHOR_VOCABULARY = [
        'CO1', 'CO2', 'CO3', 'CYTB', 'ATP6', 'ATP8',
        'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6',
        'RNR1', 'RNR2',
        'TA', 'TC', 'TD', 'TE', 'TF', 'TG', 'TH', 'TI', 'TK', 'TL1', 'TL2',
        'TM', 'TN', 'TP', 'TQ', 'TR', 'TS1', 'TS2', 'TT', 'TV', 'TW', 'TY',
    ] as Set

    static boolean isValidAnchor(anchor) {
        anchor ? (anchor.toString().trim() in ANCHOR_VOCABULARY) : false
    }

    /**
     * The gene this sample's published mitogenome is re-origined to, as a MITOS feature
     * key (CO1, RNR2, TM, TF, ...). Resolved order -> class -> group -> default.
     *
     * NEVER null: an unmapped taxon takes the table's default_anchor. This is
     * deliberately the cox1PanelGroup() trade-off rather than the seedDbGroup() one.
     * A wrong anchor rotates a circle, which is cosmetic, reversible, and softened
     * further by mitos_to_emma falling back when the gene is not annotated; no anchor
     * at all would be a hard failure. The cost of guessing here is nothing, the cost
     * of guessing a seed is a wasted reseed.
     *
     * The group step matters for a class the reference database does not cover: Calcarea
     * has no records at all, so it has no `classes` row, and without the group step a
     * calcarean sponge would take the cox1 default instead of the poriferan rrnL its
     * phylum uses 90% of the time.
     */
    static String originAnchor(taxOrder, taxClass) {
        if (!originAnchorsLoaded()) {
            throw new IllegalStateException(
                "InvertTaxonGroups.originAnchor() called before " +
                "assets/taxonomy/mito_origin_anchors.json was parsed -- call " +
                "loadOriginAnchors(file(\"\${projectDir}/assets/taxonomy/mito_origin_anchors.json\")) first")
        }
        def o = norm(taxOrder)
        def c = norm(taxClass)
        if (o && !(o in UNRESOLVED_TAXON) && originAnchorOrders.containsKey(o)) {
            return originAnchorOrders[o]
        }
        if (c && !(c in UNRESOLVED_TAXON) && originAnchorClasses.containsKey(c)) {
            return originAnchorClasses[c]
        }
        def group = seedDbGroup(c)
        if (group && originAnchorGroups.containsKey(group)) {
            return originAnchorGroups[group]
        }
        return originAnchorDefault
    }

    // Values that mean "the taxonomy did not resolve". Mirrors _UNRESOLVED in
    // bin/mito_gene_order.py and isUnresolvedTaxon in prepare_samplesheet, so the
    // Groovy and Python rank lookups behave identically on messy taxonomy.
    static final Set<String> UNRESOLVED_TAXON = ['', 'unknown', 'na', 'none', 'dropped'] as Set

    // Which curated seed database (assets/refdb/<group>/) GETORGANELLE_RESEED should
    // seed a failed invertebrate first pass from. The group name IS the directory and
    // file prefix, and must match a key of GROUPS in bin/build_invert_reference_db.py.
    //
    // Returns null for a class with no curated database, and the caller must then skip
    // the reseed and keep the first-pass assembly. There is deliberately no catch-all:
    // until now every invertebrate was reseeded from the Anthozoa database, so a
    // mollusc or a sea star was re-run against a seed far too divergent to assemble
    // from -- a guaranteed second failure that still costs the full GetOrganelle
    // walltime. Not reseeding is strictly better than reseeding from the wrong phylum.
    //
    // Cnidaria maps to 'anthozoa' because that is the only cnidarian database built so
    // far; a Hydrozoa/Scyphozoa sample is seeded from anthozoans, which is within-phylum
    // and the closest available, not cross-phylum.
    static String seedDbGroup(taxClass) {
        def c = norm(taxClass)
        if (c in CNIDARIA_CLASSES) return 'anthozoa'
        if (c in PORIFERA_CLASSES) return 'porifera'
        if (c in MOLLUSCA_CLASSES) return 'mollusca'
        if (c in ARTHROPODA_CLASSES) return 'arthropoda'
        if (c in ECHINODERMATA_CLASSES) return 'echinodermata'
        if (c in CTENOPHORA_CLASSES) return 'ctenophora'
        if (c in TUNICATA_CLASSES) return 'tunicata'
        if (c in ANNELIDA_CLASSES) return 'annelida'
        return null
    }
}
