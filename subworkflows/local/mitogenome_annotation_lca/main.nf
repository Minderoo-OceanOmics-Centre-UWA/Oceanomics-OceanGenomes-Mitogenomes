/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mitogenome modules
include { DOWNLOAD_BLAST_DB      } from '../../../modules/local/download_blast_db'
include { DOWNLOAD_TAXONKIT_DB   } from '../../../modules/local/download_taxonkit_db'
include { EMMA                   } from '../../../modules/local/EMMA'
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
    sql_config // params.sql_config
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Download taxonomy database
    //

    DOWNLOAD_BLAST_DB(Channel.value("taxdb"))
    
    // 
    // MODULE: Download taxonkit database
    //

    DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))
    
    //
    // Check if mt_assembly_prefix already exists in meta if this subworkflow is running after the mitogenome assmebly subworkflow.
    // Otherwise extract the assembly name from the fasta file and embeds it into the meta map.
    //

    fasta_with_mt_assembly_prefix = mito_assembly
    .map { meta, fasta ->
        def meta_ext = meta.containsKey('mt_assembly_prefix') 
            ? meta 
            : meta + [ mt_assembly_prefix: fasta.baseName ]
        [meta_ext, fasta]
    }
    
    //
    // MODULE: Run the program EMMA to do the annotation of the mitogenome
    //
               
    EMMA (
        fasta_with_mt_assembly_prefix // tuple val(meta), path(fasta), val(assembly_prefix)
    )

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

    combined_sequences = EMMA.out.co1_sequences
        .map { meta, file -> 
            def annotation_name = getAnnotationName(file.name)
            [meta, file, 'CO1', annotation_name] 
        }
        .mix(
            EMMA.out.s12_sequences.map { meta, file -> 
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '12s', annotation_name] 
            },
            EMMA.out.s16_sequences.map { meta, file -> 
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '16s', annotation_name] 
            }
        )
    
    //
    // MODULE: Using CO1,12s and 16s run BLAST and filter the results to provide matches for the calculation of the LCA
    //

    BLAST_BLASTN (
        combined_sequences, // tuple val(meta), path(fasta), val(gene_type), val(annotation_name)
        curated_blast_db, // params.curated_blast_db
        DOWNLOAD_BLAST_DB.out.db_files // path(db)
    )


    //
    // MODULE: Calculate the Lowest Common Ancestor (LCA) from the filtered BLAST results
    //

    LCA (
        BLAST_BLASTN.out.filtered,
        // valid_blast_results, // tuple val(meta), path(blast_filtered), val(gene_type), val(annotation_name)
        DOWNLOAD_TAXONKIT_DB.out.db_files // path(db)
    )


    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    // Need to update this section to include everything
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect{it[1]})
    ch_versions = ch_versions.mix(EMMA.out.versions.first())
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    ch_versions = ch_versions.mix(LCA.out.versions.first())



    //
    // Collate and save software versions
    //

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oceangenomes_draftgenomes_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // Emit outputs
    //

    emit:
    multiqc_files           = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    annotation_results      = EMMA.out.results
    blast_filtered_results  = BLAST_BLASTN.out.validation
    lca_results             = LCA.out.lca
    versions                = ch_collated_versions              // channel: [ path(versions.yml) ]
}
