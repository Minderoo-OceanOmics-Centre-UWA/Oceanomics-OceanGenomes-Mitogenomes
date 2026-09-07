#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Fallback genetic code only. The per-sample code is resolved from taxonomic class
// below; this is what a sample whose class has no confirmed code falls back to.
params.translation_table = params.translation_table ?: 2
// Uploads are on by default, matching the main pipeline; --skip_upload_results true
// turns this entrypoint back into the read-only QC pass it used to be.
params.skip_upload_results = params.skip_upload_results == null ? false : params.skip_upload_results

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Standalone QC-only workflow
    - Input: precomputed annotation files (*.fa/*.fasta/*.gff/*.tbl/*.gb)
    - Species: queried from SQL via VALIDATED_SPECIES_QUERY (lca_validation.validated_species_name)
    - Genetic code: resolved per sample from the taxonomic class that query also returns,
      via MitoGeneticCode.forClass(); --translation_table is the fallback for a class with
      no confirmed code (a warning names each such sample).
    - Action: run MITOGENOME_QC, then push that stage's own results to SQL via
      UPLOAD_ENA_RESULTS (ena_validation_attempts + lca_validation.validator_2).
      Only the QC stage's uploads run here -- there is no assembly, annotation or
      LCA on this path, so nothing writes mitogenome_data, blast_filtered_lca or lca.
      Disable the uploads with --skip_upload_results true.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MITOGENOME_QC } from './subworkflows/local/mitogenome_qc/main'
include { UPLOAD_ENA_RESULTS } from './subworkflows/local/upload_results_mito/main'
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
        error "Please provide --annotation_files with a glob for annotation files, e.g. --annotation_files '/path/to/mitogenomes/*/*/annotation/*.{fa,fasta,gff,tbl,gb}'"
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

            // Keep the annotation stem as the channel identity, while retaining
            // the four-field assembly prefix expected by source-modifier naming.
            def annotation_prefix = stem
            def mt_assembly_prefix = parts[0..3].join('.')
            def meta = [
                id                : parts[0],
                sequencing_type   : parts[1],
                date              : parts[2],
                code              : parts[3],
                annotation        : parts[4],
                annotation_prefix : annotation_prefix,
                mt_assembly_prefix: mt_assembly_prefix,
                // Placeholder only. The real per-sample code is resolved below from
                // the taxonomic class VALIDATED_SPECIES_QUERY returns, and overwrites
                // this before MITOGENOME_QC ever sees the meta. It is set here so that
                // anything reading the meta between here and that point (and any
                // future consumer of ch_annotations_grouped) still finds a valid code
                // rather than null.
                genetic_code      : (params.translation_table ?: 2) as int
            ]
            [ annotation_prefix, meta, file ]
        }
        .groupTuple(by: 0)
        .map { _assembly_prefix, metas, files ->
            def annotation_files = files.flatten()
            def meta = GffCircularity.annotate(metas[0], annotation_files)
            tuple(meta, annotation_files)
        }

    // Query validated species names from lca_validation.
    VALIDATED_SPECIES_QUERY(
        ch_annotations_grouped.map { meta, _files -> meta },
        sql_config_file
    )

    ch_species = VALIDATED_SPECIES_QUERY.out.species
        .map { meta, species ->
            def species_name = species ? species.toString().trim() : 'unknown'
            [ meta.annotation_prefix, species_name ?: 'unknown' ]
        }

    // The taxonomic class behind that species, for the mitochondrial genetic code.
    //
    // This entrypoint has no samplesheet, so it cannot take meta.genetic_code the way
    // the main pipeline does -- it used to assume the run-level --translation_table for
    // every sample, i.e. the vertebrate code 2 unless the operator remembered to pass
    // otherwise. That is not a QC-only concern: meta.genetic_code becomes the mgcode in
    // GEN_FILES_TABLE2ASN and the table in TRANSLATE_GENES and FORMAT_FILES, so
    // re-QCing a coral assembly here silently rewrote its annotation to code 2 and
    // submitted it that way, having been annotated under code 4 by the run that
    // produced it.
    //
    // The class comes from the same SQL round trip as the species name, so this costs
    // no extra process, and MitoGeneticCode.forClass() is the same lookup
    // prepare_samplesheet uses -- the two entrypoints cannot disagree on a class.
    ch_tax_class = VALIDATED_SPECIES_QUERY.out.tax_class
        .map { meta, tax_class_file ->
            def tax_class = tax_class_file.text?.trim() ?: ''
            [ meta.annotation_prefix, tax_class ]
        }

    // Build the tuple shape required by MITOGENOME_QC.
    ch_qc_input = ch_annotations_grouped
        .map { meta, files -> [ meta.annotation_prefix, meta, files ] }
        .join(ch_species, by: 0)
        .join(ch_tax_class, by: 0)
        .map { _annotation_prefix, meta, files, species_name, tax_class ->
            def default_code = (params.translation_table ?: 2) as int
            def mapped_code = MitoGeneticCode.forClass(tax_class)
            // Unlike prepare_samplesheet, an unmapped class does not abort here. That
            // check exists to stop a wrong table being baked into an annotation that
            // is about to be made; these annotations already exist, and refusing to
            // re-QC them helps nobody. Warn instead, naming the sample and the class,
            // so an operator can re-run with --translation_table and know which
            // samples needed it.
            if (mapped_code == null) {
                log.warn "QC-only: ${meta.annotation_prefix}: no confirmed mitochondrial " +
                         "genetic code for class '${tax_class ?: 'unresolved'}' -- using " +
                         "--translation_table ${default_code}. Pass --translation_table " +
                         "explicitly if that is wrong for this sample."
            }
            def qc_meta = meta + [ genetic_code: mapped_code ?: default_code ]
            tuple(qc_meta, species_name, true, qc_meta.circular as boolean, files)
        }

    MITOGENOME_QC(
        ch_qc_input
    )

    // Push the QC stage's own results to SQL. The scope is deliberately the QC
    // stage only: this entrypoint runs no assembly, no annotation and no LCA, so
    // it pushes no mitogenome_data, blast_filtered_lca or lca rows. That also
    // keeps SPECIES_VALIDATION out of this path, which matters -- under
    // --force_db_overwrite that module overwrites lca_validation.validated_species_name
    // and validator, and this entrypoint exists precisely for samples whose species
    // was validated by hand and cannot be re-derived from BLAST.
    //
    // prior_upload_status_files is empty here for the same reason: the five
    // pre-QC pushes never ran, so the upload report is built from the ENA
    // validation and validator_2 pushes alone.
    if (!params.skip_upload_results) {
        UPLOAD_ENA_RESULTS(
            MITOGENOME_QC.out.ena_validation_records,
            Channel.empty(),
            sql_config_file
        )
    }
}

workflow {
    QC_ONLY_FROM_ANNOTATIONS()
}
