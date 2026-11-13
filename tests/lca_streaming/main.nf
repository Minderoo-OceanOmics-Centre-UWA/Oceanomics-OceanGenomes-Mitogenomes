// Stub harness for the per-sample LCA result grouping in
// subworkflows/local/upload_results_mito.
//
// The grouping must release each sample as soon as its OWN regions are done,
// rather than waiting for every LCA task in the run. It must also cope with a
// sample producing 0, 1, 2 or 3 regions. Both properties are exercised here
// without needing the SQL processes the real subworkflow runs.

nextflow.enable.dsl = 2

// Import the real implementations rather than copying them, so this test cannot
// silently pass against a stale duplicate of the logic it is meant to protect.
include { countAnnotatedRegions     } from '../../subworkflows/local/mitogenome_annotation_lca/main.nf'
include { groupResultsByRegionCount } from '../../subworkflows/local/upload_results_mito/main.nf'

// Stands in for LCA: one task per region, with a per-sample delay so the test
// can prove a fast sample is not held behind a slow one.
process REGION_LCA {
    input:
    tuple val(meta), val(region), val(delay_s)

    output:
    tuple val(meta), path("lca.${region}.${meta.id}.tsv")

    script:
    """
    sleep ${delay_s}
    printf 'species_in_LCA\\n' > lca.${region}.${meta.id}.tsv
    """
}

// Stands in for SPECIES_VALIDATION: records when the sample was released.
process CONSUME_GROUP {
    input:
    tuple val(meta), path(lca_files)

    output:
    tuple val(meta), path("released.${meta.id}.txt")

    script:
    """
    printf '%s\\t%s\\t%s\\n' "${meta.id}" "\$(ls -1 ${lca_files} | wc -l)" "\$(date +%s%N)" \\
        > released.${meta.id}.txt
    """
}

workflow LCA_STREAMING {

    take:
    annotation_bundles // [ meta, [annotation files incl. cds/] ]
    region_work        // [ meta, region, delay_s ]

    main:
    region_counts = annotation_bundles
        .map { meta, files -> [ meta, countAnnotatedRegions(files) ] }

    REGION_LCA(region_work)

    grouped = groupResultsByRegionCount(REGION_LCA.out, region_counts.filter { _m, n -> n > 0 })

    // Zero-region samples never reach LCA, so they are injected directly rather
    // than joined -- mirroring ch_zero_region_blast_lca in the real subworkflow.
    zero_region = region_counts
        .filter { _meta, n_regions -> n_regions == 0 }
        .map { meta, _n -> [ meta, [ file("${projectDir}/assets/placeholders/empty_lca.tsv", checkIfExists: true) ] ] }

    CONSUME_GROUP(grouped.mix(zero_region))

    emit:
    counts   = region_counts
    grouped  = grouped
    released = CONSUME_GROUP.out
}
