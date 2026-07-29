nextflow.enable.dsl = 2

workflow GFF_CIRCULARITY {
    take:
    annotation_groups

    main:
    qc_input = annotation_groups.map { meta, files, species_name ->
        def annotated_meta = GffCircularity.annotate(meta, files)
        tuple(annotated_meta, species_name, true, annotated_meta.circular as boolean, files)
    }

    emit:
    qc_input
}

// Lightweight entry point for exercising the parser without nf-test.
workflow {
    if (!params.gffs) {
        error "Provide --gffs as a comma-separated list of GFF paths."
    }

    annotation_groups = Channel.fromList(params.gffs.toString().split(',') as List)
        .map { gff_path ->
            def gff = file(gff_path, checkIfExists: true)
            def prefix = gff.baseName
            tuple([id: prefix, mt_assembly_prefix: prefix], [gff], 'Test species')
        }

    GFF_CIRCULARITY(annotation_groups)
    GFF_CIRCULARITY.out.qc_input.view { meta, species, proceed, circular, files ->
        "QC_INPUT\t${meta.mt_assembly_prefix}\t${species}\t${proceed}\t${circular}\t${files.size()}"
    }
}
