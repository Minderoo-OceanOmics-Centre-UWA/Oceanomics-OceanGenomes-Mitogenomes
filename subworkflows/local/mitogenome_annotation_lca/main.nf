/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mitogenome modules
include { DOWNLOAD_BLAST_DB      } from '../../../modules/local/download_blast_db'
include { DOWNLOAD_TAXONKIT_DB   } from '../../../modules/local/download_taxonkit_db'
include { EMMA                   } from '../../../modules/local/EMMA'
include { MITOS2                 } from '../../../modules/local/mitos2'
include { BLAST_BLASTN           } from '../../../modules/nf-core/blast/blastn'
include { LCA                    } from '../../../modules/local/LCA'

// Helper functions
include { softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'


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

    MITOS2 (
        ch_annot_branched.invert.map { meta, fasta ->
            if (!mitos_refdb) {
                error "Sample ${meta.id} is marked invertebrates=true, but --mitos_refdb was not provided"
            }
            [meta, fasta]
        },
        ch_mitos_refdb
    )

    // Merge the two annotators' outputs so everything downstream is unchanged.
    ch_annot_co1     = EMMA.out.co1_sequences.mix(MITOS2.out.co1_sequences)
    ch_annot_s12     = EMMA.out.s12_sequences.mix(MITOS2.out.s12_sequences)
    ch_annot_s16     = EMMA.out.s16_sequences.mix(MITOS2.out.s16_sequences)
    ch_annot_results = EMMA.out.results.mix(MITOS2.out.results)
    ch_annot_params  = EMMA.out.tool_params.mix(MITOS2.out.tool_params)
    ch_annot_versions = EMMA.out.versions.mix(MITOS2.out.versions)

    //
    // Function to extract annotation name
    //

    def getAnnotationName = { filename ->
        def name = filename.toString().replaceAll(/\.fa$/, '')
        def parts = name.split('\\.', 2)
        return parts.size() > 1 ? parts[1] : name
    }

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



    //
    // Emit outputs
    //

    emit:
    multiqc_files           = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    annotation_results      = ch_annot_results
    blast_filtered_results  = BLAST_BLASTN.out.validation
    lca_results             = LCA.out.lca
    lca_raw_results         = LCA.out.lca_raw
    versions                = ch_versions              // channel: [ path(versions.yml) ]
}
