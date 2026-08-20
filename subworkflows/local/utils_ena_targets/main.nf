/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RESOLVE THE ENA STUDY FOR A CANDIDATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The study is created by the upstream pipeline and handed to this one as a
    single accession, so every candidate in a run submits to the same study.

    There is deliberately no default: an unset study would silently validate
    candidates against whatever the run happened to fall back to, so a missing
    or malformed accession stops the run instead.

    Locus-tag prefixes are registered against the study, but this pipeline does
    not resolve or apply them: the downstream submission pipeline owns locus
    tags end to end.

    These are functions in a module rather than a class in lib/ so that they
    resolve from any entry script, including the stub workflows under tests/.
*/

/*
    Startup guard: the run's study must be set and must look like an ENA study
    accession before any candidate is packaged.
*/
def validateEnaStudy(Map params) {
    def study = params?.ena_study?.toString()?.trim()
    if (!study) {
        throw new IllegalArgumentException(
            "--ena_study is not set, so candidates have no ENA study."
        )
    }
    if (!(study ==~ /PRJEB[0-9]+/)) {
        throw new IllegalArgumentException(
            "--ena_study must be an ENA study accession (PRJEB...), got '${study}'."
        )
    }
    return study
}

/* Add the run's study to a candidate's meta map. */
def enaStudyAnnotate(Map params, Map meta) {
    return meta + [ena_study: validateEnaStudy(params)]
}
