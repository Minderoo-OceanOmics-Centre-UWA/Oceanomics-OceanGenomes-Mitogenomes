/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mitogenome modules
include { DOWNLOAD_BLAST_DB      } from '../../../modules/local/download_blast_db'
include { DOWNLOAD_TAXONKIT_DB   } from '../../../modules/local/download_taxonkit_db'
include { EMMA                   } from '../../../modules/local/EMMA'
include { ROTATE_ORIGIN          } from '../../../modules/local/rotate_origin'
include { MITOS2                 } from '../../../modules/local/mitos2'
include { ANNOTATION_QC_GATE     } from '../../../modules/local/annotation_qc_gate'
include { CORAL_ANNOTATION_FIX   } from '../../../modules/local/coral_annotation_fix'
include { REFERENCE_RELEVANCE    } from '../../../modules/local/reference_relevance'
include { SELECT_CORAL_REFERENCE } from '../../../modules/local/select_coral_reference'
include { BLAST_BLASTN           } from '../../../modules/nf-core/blast/blastn'
include { LCA                    } from '../../../modules/local/LCA'

// Helper functions
include { softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'

// Extract the annotation name from a per-gene FASTA filename: strip the trailing
// .fa and return everything after the first dot (e.g. CO1.<prefix>.fa -> <prefix>).
// Defined at file scope so it resolves inside the .map closures below (a closure
// assigned to a local `def` is not in scope for nested operator closures).
def getAnnotationName(filename) {
    def name = filename.toString().replaceAll(/\.fa$/, '')
    def parts = name.split('\\.', 2)
    return parts.size() > 1 ? parts[1] : name
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME ANNOTATION AND LCA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ANNOTATION {

    take:
    mito_assembly //  tuple val(meta), path(fasta)
    curated_blast_db // params.curated_blast_db
    nt_blast_db // params.nt_blast_db
    mitos_refdb // params.mitos_refdb (MITOS2 RefSeq reference data dir, for invertebrates)
    reference_gb // tuple val(meta), path(reference.gb) - per-sample, partial (assembly stage)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Download taxonomy database
    //

    DOWNLOAD_BLAST_DB(Channel.value("taxdb"))
    
    ch_blast_db = DOWNLOAD_BLAST_DB.out.db_files

    // 
    // MODULE: Download taxonkit database
    //

    DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))
    
    //
    // Extract the assembly name from the fasta file and embeds it into the meta map.
    // Overwrites mt_assembly_prefix if it already exists in meta to allow for any concatgenated sequences created in sanitse fasta module.
    // The concatinated fasta is to allow for species validation from multiple contig/scaffold assemblies.
    //

    fasta_with_mt_assembly_prefix = mito_assembly
    .map { meta, fasta ->
        def meta_ext = meta + [ mt_assembly_prefix: fasta.baseName ]
        [meta_ext, fasta]
    }

    //
    // MODULE: Reference relevance check.
    // The reference is resolved from the sample's species label, so a wrong/coarse
    // label yields a wrong-family reference that silently degrades seeding and the
    // coral annotation fix. BLAST the resolved reference against the assembly and
    // record a PASS/MISMATCH review flag for every sample that has one. Label- and
    // taxonomy-DB-free; always exits 0.
    //
    REFERENCE_RELEVANCE (
        fasta_with_mt_assembly_prefix.join(reference_gb)
            .map { meta, fasta, ref -> [meta, fasta, ref] }
    )
    ch_reference_relevance = REFERENCE_RELEVANCE.out.flag

    //
    // MODULE: Annotate the mitogenome.
    // Route by meta.invertebrates: vertebrates -> EMMA, invertebrates -> MITOS2.
    // EMMA does not annotate inverts correctly, so coral/invert samples use
    // MITOS2 instead. MITOS2 emits the same output channels as EMMA so the rest
    // of this subworkflow is annotator-agnostic.
    //

    ch_annot_branched = fasta_with_mt_assembly_prefix.branch { meta, _fasta ->
        invert: meta.invertebrates
        vert:   true
    }

    // MITOS2 needs its RefSeq reference data dir; fail clearly if an invert
    // sample is present but --mitos_refdb was not provided.
    ch_mitos_refdb = mitos_refdb
        ? Channel.fromPath(mitos_refdb, checkIfExists: true).first()
        : Channel.value([])

    EMMA (
        ch_annot_branched.vert // tuple val(meta), path(fasta)
    )

    // Re-origin invert (coral) assemblies to the cox1 start before MITOS2 so the
    // anthozoan nad5 group I intron (and cox3) no longer straddle GetOrganelle's
    // linearisation point. This mirrors EMMA's `--rotate MT-TF` for vertebrates;
    // corals lack tRNA-Phe, so cox1 is the anchor. The original (un-rotated)
    // assembly is published separately by the assembly stage and left untouched.
    ch_cox1_ref = Channel.fromPath("${projectDir}/assets/cox1_anthozoa.faa", checkIfExists: true).first()

    ROTATE_ORIGIN (
        ch_annot_branched.invert.map { meta, fasta ->
            if (!mitos_refdb) {
                error "Sample ${meta.id} is marked invertebrates=true, but --mitos_refdb was not provided"
            }
            [meta, fasta]
        },
        ch_cox1_ref
    )

    MITOS2 (
        ROTATE_ORIGIN.out.fasta,
        ch_mitos_refdb
    )

    //
    // Anthozoan annotation QC gate + reference-based fixer.
    // MITOS2 annotates coral PCGs + 12S correctly but routinely drops the
    // divergent 16S and one exon of the intron-split nad5. The gate flags each
    // invert annotation as FIX (deficient) or PASS, so only broken corals are
    // re-annotated; correctly annotated batch-mates pass through MITOS2 untouched.
    //
    ANNOTATION_QC_GATE ( MITOS2.out.gff_proteins )

    ch_gate = ANNOTATION_QC_GATE.out.decision
        .map { meta, dfile -> [meta, dfile.text.split('\t')[0].trim()] }   // [meta, FIX|PASS]

    // PASS corals: keep MITOS2 output unchanged.
    ch_mitos_co1_pass     = MITOS2.out.co1_sequences.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }
    ch_mitos_s12_pass     = MITOS2.out.s12_sequences.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }
    ch_mitos_s16_pass     = MITOS2.out.s16_sequences.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }
    ch_mitos_results_pass = MITOS2.out.results.join(ch_gate).filter { it[-1] == 'PASS' }.map { meta, f, _d -> [meta, f] }

    // FIX corals: fixer base = the cox1-rotated genome MITOS2 annotated + the raw BED.
    ch_fix_base = ROTATE_ORIGIN.out.fasta
        .join(MITOS2.out.bed)
        .join(ch_gate)
        .filter { it[-1] == 'FIX' }
        .map { meta, genome, bed, _d -> [meta, genome, bed] }

    // Reference resolution -- label-free. Pick the annotation reference from the
    // curated Anthozoa mitogenome DB by SEQUENCE similarity to the assembly, so a
    // wrong/coarse species label can no longer hand a wrong-family reference to
    // the fixer. The bundled curated anthozoan reference is the only fallback (used
    // just for the rare sample no DB record aligns to).
    ch_coral_db = Channel.fromPath("${projectDir}/assets/coral_mito_refdb.gb", checkIfExists: true).first()

    SELECT_CORAL_REFERENCE ( ch_fix_base.map { meta, genome, _bed -> [meta, genome] }, ch_coral_db )
    ch_selected_ref = SELECT_CORAL_REFERENCE.out.reference   // [meta, reference.gb]

    ch_fix_selected = ch_fix_base.join(ch_selected_ref)
        .map { meta, genome, bed, ref -> [meta, genome, bed, ref] }

    // Fallback: a FIX sample for which the selector emitted no reference (no DB
    // record aligned) gets the bundled curated anthozoan reference.
    ch_sel_keys = ch_selected_ref.map { meta, _ref -> [meta, true] }
    ch_curated_ref = Channel.fromPath("${projectDir}/assets/anthozoa_reference.gb", checkIfExists: true).first()
    ch_fix_fallback = ch_fix_base.join(ch_sel_keys, remainder: true)
        .filter { it[0] != null && it[1] != null && it[-1] == null }   // FIX sample, selector emitted nothing
        .map { meta, genome, bed, _flag -> [meta, genome, bed] }
        .combine(ch_curated_ref)
        .map { meta, genome, bed, ref -> [meta, genome, bed, ref] }

    ch_coral_fix_input = ch_fix_selected.mix(ch_fix_fallback)

    CORAL_ANNOTATION_FIX ( ch_coral_fix_input )

    // Merge annotators: verts (EMMA) + PASS corals (MITOS2) + FIX corals (fixer).
    ch_annot_co1     = EMMA.out.co1_sequences.mix(ch_mitos_co1_pass, CORAL_ANNOTATION_FIX.out.co1_sequences)
    ch_annot_s12     = EMMA.out.s12_sequences.mix(ch_mitos_s12_pass, CORAL_ANNOTATION_FIX.out.s12_sequences)
    ch_annot_s16     = EMMA.out.s16_sequences.mix(ch_mitos_s16_pass, CORAL_ANNOTATION_FIX.out.s16_sequences)
    ch_annot_results = EMMA.out.results.mix(ch_mitos_results_pass, CORAL_ANNOTATION_FIX.out.results)
    ch_annot_params  = EMMA.out.tool_params.mix(ROTATE_ORIGIN.out.tool_params, MITOS2.out.tool_params, CORAL_ANNOTATION_FIX.out.tool_params)
    ch_annot_versions = EMMA.out.versions.mix(ROTATE_ORIGIN.out.versions, MITOS2.out.versions,
                                              ANNOTATION_QC_GATE.out.versions, CORAL_ANNOTATION_FIX.out.versions,
                                              SELECT_CORAL_REFERENCE.out.versions)

    //
    // Use mix() to process CO1, 12s and 16s sequences through blast
    //

    combined_sequences = ch_annot_co1
        .map { meta, file ->
            def annotation_name = getAnnotationName(file.name)
            [meta, file, 'CO1', annotation_name]
        }
        .mix(
            ch_annot_s12.map { meta, file ->
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '12s', annotation_name]
            },
            ch_annot_s16.map { meta, file ->
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '16s', annotation_name]
            }
        )
    
    //
    // Select the BLAST database per sample (curated for vertebrates, nt for invertebrates)
    //

    combined_sequences_with_db = combined_sequences
        .map { meta, file, gene_type, annotation_name ->
            if (meta.invertebrates && !nt_blast_db) {
                error "Sample ${meta.id} is marked invertebrates=true, but --nt_blast_db was not provided"
            }
            def blast_db = (meta.invertebrates ? nt_blast_db : curated_blast_db)
            [meta, file, gene_type, annotation_name, blast_db]
        }
    
    //
    // MODULE: Using CO1,12s and 16s run BLAST and filter the results to provide matches for the calculation of the LCA
    //

    BLAST_BLASTN (
        combined_sequences_with_db, // tuple val(meta), path(fasta), val(gene_type), val(annotation_name), val(blast_db)
        ch_blast_db // path(db)
    )


    //
    // MODULE: Calculate the Lowest Common Ancestor (LCA) from the filtered BLAST results
    //
    ch_worms = channel.fromPath("https://raw.githubusercontent.com/Minderoo-OceanOmics-Centre-UWA/LCA_With_Fishbase/main/data/worms_species.txt.gz")

    LCA (
        BLAST_BLASTN.out.filtered,
        ch_worms.first()
        // valid_blast_results, // tuple val(meta), path(blast_filtered), val(gene_type), val(annotation_name)
        // DOWNLOAD_TAXONKIT_DB.out.db_files // path(db)
    )


    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    // Need to update this section to include everything
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(ch_annot_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(LCA.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(ch_annot_versions.first())
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    ch_versions = ch_versions.mix(LCA.out.versions.first())
    ch_versions = ch_versions.mix(REFERENCE_RELEVANCE.out.versions.first())



    //
    // Emit outputs
    //

    emit:
    multiqc_files           = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    annotation_results      = ch_annot_results
    blast_filtered_results  = BLAST_BLASTN.out.validation
    lca_results             = LCA.out.lca
    lca_raw_results         = LCA.out.lca_raw
    reference_relevance     = ch_reference_relevance   // channel: [ meta, path(reference_relevance.txt) ] (PASS|MISMATCH|UNKNOWN)
    versions                = ch_versions              // channel: [ path(versions.yml) ]
}
