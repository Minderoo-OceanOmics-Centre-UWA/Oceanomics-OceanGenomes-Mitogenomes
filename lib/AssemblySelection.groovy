/**
 * Which of two assemblies of the same sample should be the one that goes forward.
 *
 * The rules used to live only in the GetOrganelle subworkflow, where they decided
 * reseed vs first pass. The precomputed ("skip assembly") path needs the same
 * decision -- without it every published variant is treated as its own assembly run
 * and they ALL reach annotation, which is how NOVA_260909_LA_mitogenomes-4 annotated
 * 110 assemblies for 80 samples instead of ~74. Sharing one implementation is
 * deliberate: the duplicated-logic failure mode is exactly how bin/mito_run_snapshot.py
 * came to disagree with the pipeline about broken CDS.
 *
 * Rules, in order:
 *   1. an empty candidate never wins;
 *   2. an empty incumbent always loses;
 *   3. circular beats non-circular -- the only outright success available here;
 *   4. a molecule the circularity check called LENGTH-ANOMALOUS loses to one it did not.
 *      This sits ahead of size on purpose: the size rule rewards the bigger molecule, and a
 *      length anomaly is the finding that the molecule is bigger than it should be. OG2936's
 *      welded first pass is 31,154 bp against a 20,312 bp reference (ratio 1.534, "review for
 *      duplication / NUMT") and beat a clean 19,791 bp single contig purely on that size. Only
 *      an explicit anomaly demotes: unknown is not treated as anomalous, so this can never
 *      promote a molecule on the absence of a check;
 *   5. a large size difference decides on size, because a variant that recovered
 *      several times more sequence in more pieces is not a regression;
 *   6. within the same size class, FEWER SPACER Ns wins. A reference-guided join
 *      manufactures one record by welding contigs with 100bp N runs, and any gene
 *      crossing a weld is corrupted: OG2998's CYTB reads ATT..TAA in the plain
 *      concatenation and ACT..TCT across the join. Rule 6 would otherwise reward
 *      the join for the single record the welding bought it;
 *   7. within the same size class, fewer records wins;
 *   8. otherwise keep the candidate, which is the more informed assembly.
 *
 * Rule 6 is a no-op on the live reseed-vs-first-pass decision, where both sides are
 * raw GetOrganelle output and carry no spacer Ns. It only engages once a welded molecule is in
 * the comparison. Rule 4 is likewise inert there: GETORGANELLE_CHECK runs AFTER this decision
 * on the live path, so no anomaly verdict exists yet. It is the precomputed path, which
 * reloads a previous run's verdicts, that has the evidence to act on.
 */
class AssemblySelection {

    /**
     * Measure an assembly FASTA: base pairs, record count, and N count.
     *
     * Unreadable or absent is reported as empty rather than thrown, because a
     * selection must still resolve when one candidate's file is missing.
     */
    static Map measure(fasta, Boolean circular = null, Boolean lengthAnomaly = null) {
        int bp = 0, records = 0, ns = 0
        try {
            fasta.text.readLines().each { line ->
                if (line.startsWith('>')) {
                    records++
                } else {
                    String seq = line.trim()
                    bp += seq.length()
                    ns += seq.count('N') + seq.count('n')
                }
            }
        } catch (ignored) {
            return [bp: 0, records: 0, ns: 0, circular: circular, length_anomaly: lengthAnomaly]
        }
        return [bp: bp, records: records, ns: ns, circular: circular, length_anomaly: lengthAnomaly]
    }

    /** True when `candidate` should displace `incumbent`. Both from measure(). */
    static boolean prefer(Map candidate, Map incumbent, double tolerance) {
        if (candidate == null || candidate.bp == 0) return false
        if (incumbent == null || incumbent.bp == 0) return true

        // Only an explicit true/false is evidence; unknown (null) must not be read
        // as non-circular, or an assembly with no check would lose on a fact nobody
        // established.
        if (candidate.circular != null && incumbent.circular != null
            && candidate.circular != incumbent.circular) {
            return candidate.circular as boolean
        }

        // A molecule KNOWN to be length-anomalous must not win the size rule below, which is
        // the rule that would reward it for the very excess the anomaly describes. Only an
        // explicit true demotes, so an unchecked molecule is never promoted by this.
        if (candidate.length_anomaly != incumbent.length_anomaly) {
            if (candidate.length_anomaly == true) return false
            if (incumbent.length_anomaly == true) return true
        }

        double ratio = (candidate.bp as double) / (incumbent.bp as double)
        if (ratio >= tolerance) return true
        if (ratio <= 1.0d / tolerance) return false

        if (candidate.ns != incumbent.ns) return candidate.ns < incumbent.ns
        if (candidate.records != incumbent.records) return candidate.records < incumbent.records
        return true
    }
}
