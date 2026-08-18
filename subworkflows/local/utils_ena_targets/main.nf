/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RESOLVE THE ENA STUDY FOR A CANDIDATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PRJEB110568 is an umbrella project and can never receive data, so each
    sequencing technology submits to its own ENA child study.

    Every technology with a viable mitogenome is published, so there is no
    fallback and no default: an unrecognised technology means the candidate
    would be submitted to the wrong study, and that must stop the run rather
    than be quietly resolved.

    Locus-tag prefixes are registered against these same child studies, but this
    pipeline does not resolve or apply them: the downstream submission pipeline
    owns locus tags end to end.

    These are functions in a module rather than a class in lib/ so that they
    resolve from any entry script, including the stub workflows under tests/.
*/

def enaTechnologies() {
    return ['hifi', 'hic', 'ilmn']
}

/*
    The technology of a candidate, taken from its meta or, failing that, from
    position 1 of its full SeqID (OG<n>.<tech>.<date>.<code>...).

    Deliberately strict: no trimming, no case folding.  mitogenome_data holds at
    least one 'hifi ' row with a trailing space, and a value like that is a
    data-entry fault that has to surface here, not be silently repaired into a
    submission.
*/
def enaTechOf(Map meta) {
    def tech = meta?.sequencing_type
    if (tech == null) {
        def seqid = (meta?.full_seqid ?: meta?.mt_assembly_prefix)?.toString()
        def parts = seqid ? seqid.tokenize('.') : []
        tech = parts.size() > 1 ? parts[1] : null
    }
    if (!(tech instanceof String) || !enaTechnologies().contains(tech)) {
        def label = (meta?.full_seqid ?: meta?.mt_assembly_prefix ?: meta?.id ?: 'candidate').toString()
        throw new IllegalArgumentException(
            "Cannot resolve ENA study for ${label}: sequencing technology '${tech}' is not " +
            "one of ${enaTechnologies().join(', ')}. Fix the source record rather than " +
            "defaulting, or the candidate is submitted to the wrong ENA study."
        )
    }
    return tech
}

/* Study accession for one technology. */
def enaTargetForTech(Map params, String tech) {
    if (!enaTechnologies().contains(tech)) {
        throw new IllegalArgumentException(
            "Unknown sequencing technology '${tech}'. Known technologies: ${enaTechnologies().join(', ')}."
        )
    }
    def study = params?."ena_study_${tech}"?.toString()?.trim()
    if (!study) {
        throw new IllegalArgumentException(
            "--ena_study_${tech} is not set, so ${tech} candidates have no ENA study."
        )
    }
    if (!(study ==~ /PRJEB[0-9]+/)) {
        throw new IllegalArgumentException(
            "--ena_study_${tech} must be an ENA study accession (PRJEB...), got '${study}'."
        )
    }
    return [ena_study: study]
}

/* Add the resolved study to a candidate's meta map. */
def enaTargetAnnotate(Map params, Map meta) {
    return meta + enaTargetForTech(params, enaTechOf(meta))
}

/*
    Startup guard: every technology must resolve, and two technologies must never
    share a study.  A duplicate here submits one technology's records under
    another's study.
*/
def validateEnaTargets(Map params) {
    def studies = [:]
    enaTechnologies().each { tech ->
        def target = enaTargetForTech(params, tech)
        def clash = studies[target.ena_study]
        if (clash) {
            throw new IllegalArgumentException(
                "ENA study ${target.ena_study} is set for both ${clash} and ${tech}. " +
                "Each technology needs its own child study."
            )
        }
        studies[target.ena_study] = tech
    }
}
