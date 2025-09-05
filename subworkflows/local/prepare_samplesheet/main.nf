//
// Subworkflow with functionality specific to the nf-core/oceangenomesmitogenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_SAMPLESHEET      } from '../../../modules/local/create_samplesheet'
include { HIC_DATE_QUERY          } from '../../../modules/local/hic_date_query'
include { samplesheetToList         } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_SAMPLESHEET {

    take:
    input               //  string: Path to input samplesheet
    input_dir           //  string: Path to input samplesheet
    samplesheet_prefix  // 
    
    main:

   // Handle different input types
    if (input && input_dir) {
        error "Please specify either --input (samplesheet) OR --input_dir (directory), not both"
    }
    
    if (input && input_dir) {
        error "Please specify either --input (samplesheet) or --input_dir (directory)"
    }
    
    if (input_dir) {
        // Resolve the glob pattern to actual files
        input_files_ch = Channel.fromPath(input_dir, checkIfExists: true)
            .collect()

        // Create samplesheet from directory
        CREATE_SAMPLESHEET(
            input_files_ch,
            "${samplesheet_prefix}.csv"
        )

        samplesheet_ch = CREATE_SAMPLESHEET.out.samplesheet
            .map { samplesheet_file ->
                samplesheetToList(samplesheet_file, "${projectDir}/assets/schema_input.json")
            }
            .flatMap { sample_list ->
                // Convert each sample record to the expected format
                sample_list.collect { sample_record ->
                    def meta = sample_record[0]
                    def fastq_1 = sample_record[1]
                    def fastq_2 = sample_record[2]
                    def sequencing_type = sample_record[3]
                    
                    // Add sequencing_type to meta
                    def enhanced_meta = meta + [sequencing_type: sequencing_type]
                    
                    if (!fastq_2 || fastq_2.toString() == '[]') {
                        return [ enhanced_meta.id, enhanced_meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ enhanced_meta.id, enhanced_meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
                }
            }
            .groupTuple()
            .map { samplesheet ->
                return validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }

    } else {
        samplesheet_ch = Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { item ->
                println "DEBUG ELSE 1: Item = ${item}"
                return item
            }
             .flatMap { sample_list ->
                // Convert each sample record to the expected format
                sample_list.collect { sample_record ->
                    def meta = sample_record[0]
                    def fastq_1 = sample_record[1]
                    def fastq_2 = sample_record[2]
                    def sequencing_type = sample_record[3]
                    
                    // Add sequencing_type to meta
                    def enhanced_meta = meta + [sequencing_type: sequencing_type]
                    
                    if (!fastq_2 || fastq_2.toString() == '[]') {
                        return [ enhanced_meta.id, enhanced_meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ enhanced_meta.id, enhanced_meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
                }
            }
            .groupTuple()
            .map { samplesheet ->
                return validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
    }

    // DEBUG
    samplesheet_ch.view { meta, reads -> "Sample: ${meta.id}, Type: ${meta.sequencing_type}, Meta: ${meta}" }
    
    // Branch by sequencing type
    samplesheet_ch.branch { meta, reads ->
        hifi: meta.sequencing_type == 'hifi'
            return [meta, reads]
        ilmn: meta.sequencing_type == 'ilmn'
            return [meta, reads]
        hic: meta.sequencing_type == 'hic'
            return [meta, reads]
    }.set { branched_samples }

    // Process HiFi data
    ch_hifi_with_meta = branched_samples.hifi
        .map { meta, reads ->
            // Extract date from HiFi filename (use first file if multiple)
            // Example: OG111_m84154_241004_105305_s3.hifi_reads.bc2068.filt.fastq.gz
            // Date is in format YYMMDD after the second underscore
            def filename = reads[0].name
            def parts = filename.split('_')
            def date = parts.size() >= 3 ? parts[2] : 'unknown'
            
            // Create assembly prefix: ${sample}.${sequencing_type}.${date}
            def assembly_prefix = "${meta.id}.${meta.sequencing_type}.${date}"
            
            // Enhanced meta map
            def enhanced_meta = meta + [
                date: date,
                assembly_prefix: assembly_prefix
            ]
            
            [enhanced_meta, reads]
        }

    // Process Illumina data  
    ch_ilmn_with_meta = branched_samples.ilmn
        .map { meta, reads ->
            // Extract date from Illumina filename (use first file)
            // Example: OG764.ilmn.240716.R1.fastq.gz
            // Date is the 3rd segment when split by '.'
            def filename = reads[0].name
            def parts = filename.split('\\.')
            def date = parts.size() >= 3 ? parts[2] : 'unknown'
            
            // Assembly prefix is first 3 parts separated by '.'
            // Example: OG764.ilmn.240716
            def assembly_prefix = parts.size() >= 3 ? 
                "${parts[0]}.${parts[1]}.${parts[2]}" : 
                "${meta.id}.${meta.sequencing_type}.${date}"
            
            // Enhanced meta map
            def enhanced_meta = meta + [
                date: date,
                assembly_prefix: assembly_prefix
            ]
            
            [enhanced_meta, reads]
        }


    // Query database to get hic sample sequencing date
    ch_query_meta = branched_samples.hic
        .map { meta, reads -> 
            [meta.id, meta] 
        }
    
    HIC_DATE_QUERY(
        ch_query_meta,
        params.sql_config
    )
    
    // Combine database results with original files
    ch_hic_with_files = branched_samples.hic
        .map { meta, reads -> 
            [meta, reads] 
        }
        .combine(
            HIC_DATE_QUERY.out.date.map { meta, date -> 
                [meta, date] 
            }, 
            by: 0
        )
        .map { meta, reads, date ->
            // Clean up the sample ID: remove everything after first set of numbers
            // Example: OG111M-3_HICL -> OG111
            def cleaned_id = meta.id.replaceAll(/^(OG\d+).*/, '$1')
            
            // Create assembly prefix: ${cleaned_sample}.${sequencing_type}.${date}
            def assembly_prefix = "${cleaned_id}.${meta.sequencing_type}.${date}"
            
            // Enhanced meta map with cleaned ID
            def enhanced_meta = meta + [
                id: cleaned_id,           // Update the ID to cleaned version
                original_id: meta.id,     // Keep original ID for reference
                date: date,
                assembly_prefix: assembly_prefix
            ]
            
            [enhanced_meta, reads]
        }

    // Combine all processed data
    ch_getorg = ch_ilmn_with_meta.mix(ch_hic_with_files)
    
    emit:
    samplesheet = samplesheet_ch
    getorg_input = ch_getorg
    mitohifi_input = ch_hifi_with_meta
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (sample_id, metas, fastqs) = input

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}