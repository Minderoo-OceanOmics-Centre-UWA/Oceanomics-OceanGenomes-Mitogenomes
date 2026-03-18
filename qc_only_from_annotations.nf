#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Ensure modules relying on this parameter have a stable default in standalone runs.
params.translation_table = params.translation_table ?: 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Standalone QC-only workflow
    - Input: precomputed annotation files (*.fa/*.fasta/*.gff/*.tbl/*.gb)
    - Species: queried from SQL via VALIDATED_SPECIES_QUERY (lca_validation.validated_species_name)
    - Action: run MITOGENOME_QC only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MITOGENOME_QC } from './subworkflows/local/mitogenome_qc/main'
include { VALIDATED_SPECIES_QUERY } from './modules/local/validated_species_query/main'

/*
Expected annotation filename prefix:
  <sample>.<sequencing_type>.<date>.<code>.<annotation>[...].<ext>
Example:
  OG764.ilmn.240716.getorg1770.emma100.gff
*/

workflow QC_ONLY_FROM_ANNOTATIONS {

    main:

    if (!params.annotation_files) {
        error "Please provide --annotation_files with a glob for annotation files, e.g. --annotation_files '/path/to/mitogenomes/*/*/emma/*.{fa,fasta,gff,tbl,gb}'"
    }
    if (!params.sql_config) {
        error "Please provide --sql_config for validated species lookup and source modifier generation."
    }
    if (!params.template_sbt) {
        error "Please provide --template_sbt for GEN_FILES_TABLE2ASN."
    }
    if (!params.outdir) {
        error "Please provide --outdir."
    }

    def sql_config_file = file(params.sql_config, checkIfExists: true)
    def template_sbt_file = file(params.template_sbt, checkIfExists: true)
    // Ensure downstream modules that consume params.* path inputs receive File objects.
    params.sql_config = sql_config_file
    params.template_sbt = template_sbt_file

    ch_annotation_files = Channel.fromPath(params.annotation_files, checkIfExists: true)

    // Group input annotation files by assembly prefix and reconstruct minimal metadata.
    ch_annotations_grouped = ch_annotation_files
        .map { file ->
            def stem = file.baseName
            def parts = stem.split('\\.')

            if (parts.length < 5) {
                error "Input file '${file.name}' must include '<sample>.<sequencing_type>.<date>.<code>.<annotation>'"
            }

            def mt_assembly_prefix = stem
            def meta = [
                id                : parts[0],
                sequencing_type   : parts[1],
                date              : parts[2],
                code              : parts[3],
                annotation        : parts[4],
                mt_assembly_prefix: mt_assembly_prefix
            ]
            [ mt_assembly_prefix, meta, file ]
        }
        .groupTuple(by: 0)
        .map { _assembly_prefix, metas, files ->
            tuple(metas[0], files.flatten())
        }

    // Query validated species names from lca_validation.
    VALIDATED_SPECIES_QUERY(
        ch_annotations_grouped.map { meta, _files -> meta },
        sql_config_file
    )

    ch_species = VALIDATED_SPECIES_QUERY.out.species
        .map { meta, species ->
            def species_name = species ? species.toString().trim() : 'unknown'
            [ meta.mt_assembly_prefix, species_name ?: 'unknown' ]
        }

    // Build the tuple shape required by MITOGENOME_QC.
    ch_qc_input = ch_annotations_grouped
        .map { meta, files -> [ meta.mt_assembly_prefix, meta, files ] }
        .join(ch_species, by: 0)
        .map { _prefix, meta, files, species_name ->
            tuple(meta, species_name, 'true', files)
        }

    MITOGENOME_QC(
        ch_qc_input
    )
}

workflow {
    QC_ONLY_FROM_ANNOTATIONS()
}
