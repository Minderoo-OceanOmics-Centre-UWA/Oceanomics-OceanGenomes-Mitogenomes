/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Pipeline subworkflows
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oceangenomesmitogenomes_pipeline'

// Mitogenome assembly subworkflows
include { MITOGENOME_ASSEMBLY_GETORG       } from '../subworkflows/local/mitogenome_assembly/getorganelle'
include { MITOGENOME_ASSEMBLY_MITOHIFI     } from '../subworkflows/local/mitogenome_assembly/mitohifi'
include { MITOGENOME_ANNOTATION     } from '../subworkflows/local/mitogenome_annotation_lca'
include { UPLOAD_RESULTS            } from '../subworkflows/local/upload_results_mito'
include { MITOGENOME_QC             } from '../subworkflows/local/mitogenome_qc'
include { SANITISE_FASTA           } from '../modules/local/sanitise_fasta/main'

// Using local module SANITISE_FASTA (see modules/local/sanitise_fasta)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OCEANGENOMESMITOGENOMES {

    take:
    getorg_input    // tuple for getorganelle - tuple(meta, reads)
    mitohifi_input // tuple for mitohifi - tuple(meta, reads)
    curated_blast_db // params.curated_blast_db
    sql_config // params.sql_config
    organelle_type // params.organelle_type "animal_mt"

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /* The if and else statements in this workflow are for when steps are skipped in the nextflow_run script.
        What it is doing is 'if' this processes isnt skipped then run the subworkflow and provide the standard outputs.
        
        Then 'else if' is if the process has been skipped then we will create the required input for the following
        subworkflows using the predefined file paths for the output in the nextflow.config, the "precomputed_*" file
        paths. 
            The relevant information is extracted to create the meta map from the parts of the file name. This is
            assuming that files are named with the sample id, then the type of sequencing, the date and then the other
            information in the file name, all seperated by a '.'
            The mitogenome sections assume that the files are named as above, $sample_id.$sequencing_technology.$date with
            additional information added to this with each proccess. The assembly process will add a .getorganelle${version}
            after the inital prefix and befor the file extension. Then the annotation process will add a .emma${version} to 
            the previous prefix.
        
        The 'else' statement is then just creating and empty channel if neither of the previous steps worked.
    */

    
    // 
    // SUBWORKFLOW: MITOGENOME_ASSEMBLY_GETORG
    //

    if (!params.skip_mitogenome_assembly_getorg) {
        MITOGENOME_ASSEMBLY_GETORG (
            getorg_input,
            organelle_type // params.organelle_type
        )
        ch_mitogenome_getorg_assembly_fasta = MITOGENOME_ASSEMBLY_GETORG.out.assembly_fasta
        ch_mitogenome_getorg_assembly_log = MITOGENOME_ASSEMBLY_GETORG.out.assembly_log
    } else if (params.precomputed_mitogenome_assembly_fasta_getorg) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_getorg_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_getorg, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the filename (without extension)
            def filename = file.baseName  // Gets filename without .fasta extension
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def date = parts.length > 2 ? parts[2] : null
            def sample_id = parts.length > 2 ? parts[0..3].join('.') : filename
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                mt_assembly_prefix: sample_id
            ]
            return tuple(meta, file)
        }
        ch_mitogenome_getorg_assembly_log = Channel.fromPath(params.precomputed_mitogenome_assembly_log_getorg, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the filename (without extension)
            def filename = file.baseName  // Gets filename without extension
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def date = parts.length > 2 ? parts[2] : null
            def sample_id = parts.length > 2 ? parts[0..3].join('.') : filename
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                mt_assembly_prefix: sample_id
            ]
            return tuple(meta, file)
        }
    } else {
        ch_mitogenome_getorg_assembly_fasta = Channel.empty()
        ch_mitogenome_getorg_assembly_log = Channel.empty()
    }


    // 
    // SUBWORKFLOW: MITOGENOME_ASSEMBLY_MITOHIFI
    //

    if (!params.skip_mitogenome_assembly_hifi) {
        MITOGENOME_ASSEMBLY_MITOHIFI (
            mitohifi_input,
        )
        ch_mitogenome_hifi_assembly_fasta = MITOGENOME_ASSEMBLY_MITOHIFI.out.assembly_fasta
        ch_mitogenome_hifi_assembly_log = MITOGENOME_ASSEMBLY_MITOHIFI.out.assembly_log
    } else if (params.precomputed_mitogenome_assembly_fasta_hifi) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_hifi_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta_hifi, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the filename (without extension)
            def filename = file.baseName  // Gets filename without .fasta extension
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def date = parts.length > 2 ? parts[2] : null
            def sample_id = parts.length > 2 ? parts[0..3].join('.') : filename
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                mt_assembly_prefix: sample_id
            ]
            return tuple(meta, file)
        }
        ch_mitogenome_hifi_assembly_log = Channel.fromPath(params.precomputed_mitogenome_assembly_log_hifi, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the filename (without extension)
            def filename = file.baseName  // Gets filename without extension
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def date = parts.length > 2 ? parts[2] : null
            def sample_id = parts.length > 2 ? parts[0..3].join('.') : filename
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                mt_assembly_prefix: sample_id
            ]
            return tuple(meta, file)
        }
    } else {
        ch_mitogenome_hifi_assembly_fasta = Channel.empty()
        ch_mitogenome_hifi_assembly_log = Channel.empty()
    }

    //
    // mix the hifi and getorganelle outputs to feed into the annotation
    //

    ch_annotation_input = ch_mitogenome_hifi_assembly_fasta.mix(ch_mitogenome_getorg_assembly_fasta)

    // Sanitise FASTA before annotation to avoid duplicate IDs / multi-contig issues
    SANITISE_FASTA(
        ch_annotation_input
    )
    ch_annotation_input_sanitised = SANITISE_FASTA.out

    //
    // SUBWORKFLOW: MITOGENOME_ANNOTATION
    //
    if (!params.skip_mitogenome_annotation) {
        MITOGENOME_ANNOTATION (
            ch_annotation_input_sanitised,
            curated_blast_db
        )
        ch_mitogenome_annotation_results = MITOGENOME_ANNOTATION.out.annotation_results
        ch_mitogenome_blast_results = MITOGENOME_ANNOTATION.out.blast_filtered_results
        ch_mitogenome_lca_results = MITOGENOME_ANNOTATION.out.lca_results
    } else if (params.precomputed_mitogenome_annotation_results) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_annotation_results = Channel.fromPath(params.precomputed_mitogenome_annotation_results)
        ch_mitogenome_blast_results = Channel.fromPath(params.precomputed_mitogenome_blast_results)
        ch_mitogenome_lca_results = Channel.fromPath(params.precomputed_mitogenome_lca_results)
    } else {
        ch_mitogenome_annotation_results = Channel.empty()
        ch_mitogenome_blast_results = Channel.empty()
        ch_mitogenome_lca_results = Channel.empty()
    }

    //
    // Combine outputs for data uploads
    //add in mitohifi stuff too

    ch_mitogenome_getorg_assembly_results = ch_mitogenome_getorg_assembly_fasta.join(ch_mitogenome_getorg_assembly_log, by: 0)
    ch_mitogenome_hifi_assembly_results = ch_mitogenome_hifi_assembly_fasta.join(ch_mitogenome_hifi_assembly_log, by: 0)

    ch_mitogenome_assembly_results = ch_mitogenome_getorg_assembly_results.mix(ch_mitogenome_hifi_assembly_results)

    //
    // SUBWORKFLOW: UPLOAD_RESULTS
    //
    // Conditional uploading of results to SQL and species check - only run if not skipped
    // All these processes access the OceanOmics PostgreSQL database.
    if (!params.skip_upload_results) {
        UPLOAD_RESULTS (
            ch_mitogenome_assembly_results,
            ch_mitogenome_annotation_results,
            ch_mitogenome_blast_results,
            ch_mitogenome_lca_results,
            sql_config // params.sql_config
        )
    }


    
    ch_qc_input = UPLOAD_RESULTS.out.qc_ready.join(ch_mitogenome_annotation_results, by:0)

    // If the LCA validation is correct, then run the QC to prepare for submission to GenBank
    // Need to add this into the pipeline.
    // Now that protein lengths are being added to the database it could provide a list of 
    // non submitted mitogenomes they can be grouped with to submit and then say when there is a 
    // group of similar mitogenomes they can be submitted as a batch.
    MITOGENOME_QC (
        ch_qc_input // tuple val(meta), val(species_name), val(proceed_qc true/false), path(emma/*)
    )

    //
    // Collect all MultiQC files from all subworkflows
    //

    if (!params.skip_mitogenome_assembly_getorg) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ASSEMBLY_GETORG.out.multiqc_files)}
    if (!params.skip_mitogenome_assembly_hifi) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ASSEMBLY_MITOHIFI.out.multiqc_files)}
    if (!params.skip_mitogenome_annotation) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ANNOTATION.out.multiqc_files)}
    if (!params.skip_upload_results) {ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_RESULTS.out.multiqc_files)}
    if (!params.skip_upload_results) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_QC.out.multiqc_files)}

    // 
    // Collect all versions from subworkflows
    //

    if (!params.skip_mitogenome_assembly_getorg) {
        ch_versions = ch_versions.mix(
            MITOGENOME_ASSEMBLY_GETORG.out.versions,
        )
    }
    if (!params.skip_mitogenome_assembly_hifi) {
        ch_versions = ch_versions.mix(
            MITOGENOME_ASSEMBLY_MITOHIFI.out.versions
        )
    }
    if (!params.skip_mitogenome_annotation) {
        ch_versions = ch_versions.mix(MITOGENOME_ANNOTATION.out.versions)
    }
    // Upload + QC subworkflows provide versions; guard independently
    if (!params.skip_upload_results) {
        ch_versions = ch_versions.mix(UPLOAD_RESULTS.out.versions)
    }
    // Run of MITOGENOME_QC depends on upstream evaluation; include if present
    ch_versions = ch_versions.mix(MITOGENOME_QC.out.versions)







    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oceangenomesmitogenomes_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
