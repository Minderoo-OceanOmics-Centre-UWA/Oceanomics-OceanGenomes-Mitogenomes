// Make a precomputed annotation bundle complete and consistent with its own GFF, for the
// --skip_mitogenome_annotation --blast_precomputed_annotations path. A file glob cannot
// carry the bundle's cds/ and proteins/ directories, and without them ANNOTATION_STATS
// skips the PCG ORF check, so both are rebuilt here; the CO1/12S/16S BLAST inputs come
// out of the same pass. Outputs mirror the annotators: `bundle` stands in for their
// path("annotation/*") results and co1/s12/s16 for their region emits. See
// bin/rebuild_annotation_bundle.py for why the GFF, not the bundle's cds/, is authoritative.
process REBUILD_ANNOTATION_BUNDLE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(annotation_files, stageAs: 'bundle/*')

    output:
    tuple val(meta), path("annotation/*"),        emit: bundle
    tuple val(meta), path("regions/CO1.*.fa"),    emit: co1_sequences, optional: true
    tuple val(meta), path("regions/RNR1.*.fa"),   emit: s12_sequences, optional: true
    tuple val(meta), path("regions/RNR2.*.fa"),   emit: s16_sequences, optional: true
    // Always emitted, one per assembly (0-3), so the precomputed path's region_counts
    // stays total over the annotation rows that MITOGENOME_QC joins against.
    tuple val(meta), env('N_REGIONS'),            emit: region_count
    tuple val(meta), path("regions/*.bundle_rebuild.tsv"), emit: status
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Only used when a protein has to be re-translated (stale or missing cds/).
    def gcode_arg = meta.genetic_code ? "--genetic-code ${meta.genetic_code}" : ''
    """
    rebuild_annotation_bundle.py --outdir annotation --regions-dir regions ${gcode_arg} bundle/*
    N_REGIONS=\$(find regions -maxdepth 1 -name '*.fa' | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def stem = "${meta.mt_assembly_prefix}.stub"
    """
    mkdir -p annotation/cds annotation/proteins regions
    cp -L bundle/* annotation/ 2>/dev/null || true
    for region in CO1 RNR1 RNR2; do
        printf '>${stem}\\nACGT\\n' > regions/\${region}.${stem}.fa
        cp regions/\${region}.${stem}.fa annotation/cds/
    done
    printf 'annotation\\tgene\\tkind\\tstatus\\tlength\\tcds_check\\tnotes\\n' > regions/${stem}.bundle_rebuild.tsv
    N_REGIONS=3
    : > versions.yml
    """
}
