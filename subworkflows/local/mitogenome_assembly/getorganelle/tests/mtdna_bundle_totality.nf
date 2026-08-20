// Test-only harness guarding the Plan 2f collapsed-variant mirror bundle invariant for
// the GetOrganelle and oatk arms (subworkflows/local/mitogenome_assembly/{getorganelle,
// mitohifi}/main.nf `ch_mtdna_files`).
//
// The bundle used to end in a bare groupTuple (emits nothing until the source closes = a
// whole-run barrier). It is now sized with groupKey so each sample's bundle releases as
// soon as its own files arrive. For GetOrganelle the count is EXACTLY 3 per sample
// (assembly + log + circularity evidence, all from ch_getorg_check_in) and for oatk
// EXACTLY 2 (contig + circularity evidence) -> plain groupTuple, no remainder, every sample
// streams. This harness reproduces those exact constructions and asserts: one bundle per
// sample, each of the fixed size. If a subworkflow's bundle contributor set changes, update
// it and this harness together.

workflow MTDNA_BUNDLE_TOTALITY {

    main:
    def f = file("${projectDir}/assets/NO_REFERENCE.gb", checkIfExists: true)   // arbitrary stand-in file

    // --- GetOrganelle bundle: fasta + log + evidence, keyed by lineage prefix, size 3 ---
    def gmeta = (1..5).collect { [ id: "OG${it}", mt_assembly_run_prefix: "G${it}" ] }
    getorg_fasta    = Channel.fromList( gmeta.collect { [ it, f ] } )
    getorg_log      = Channel.fromList( gmeta.collect { [ it, f ] } )
    getorg_evidence = Channel.fromList( gmeta.collect { [ it, f ] } )
    getorg_bundle = getorg_fasta
        .mix( getorg_log, getorg_evidence )
        .map { meta, x -> [ groupKey(meta.mt_assembly_run_prefix, 3), x ] }
        .groupTuple()
        .map { key, files -> [ key.getGroupTarget(), files.size() ] }

    // --- oatk bundle: contig + circularity evidence, size 2 ---
    def ometa = (1..4).collect { [ id: "OG${it}", mt_assembly_run_prefix: "K${it}" ] }
    oatk_fasta    = Channel.fromList( ometa.collect { [ it, f ] } )
    oatk_evidence = Channel.fromList( ometa.collect { [ it, f ] } )
    oatk_bundle = oatk_fasta
        .mix( oatk_evidence )
        .map { meta, x -> [ groupKey(meta.mt_assembly_run_prefix, 2), x ] }
        .groupTuple()
        .map { key, files -> [ key.getGroupTarget(), files.size() ] }

    emit:
    getorg_bundle = getorg_bundle
    oatk_bundle   = oatk_bundle
}
