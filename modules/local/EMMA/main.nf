process EMMA {
    tag "$meta.id" // tag is is being used for publish dir
    label 'process_medium'

    // No native Singularity container available, using Docker image for both engines
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/emma:1.0' :
        'tylerpeirce/emma:1.0' }"
    
    

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("emma/cds/*CO1*.fa"),  emit: co1_sequences
    tuple val(meta), path("emma/cds/*RNR1*.fa"),  emit: s12_sequences  
    tuple val(meta), path("emma/cds/*RNR2*.fa"),  emit: s16_sequences
    tuple val(meta), path("emma/*"), emit: results 
    tuple val(meta), path("07_emma.tool_params_mqcrow.html"), emit: tool_params
    path "versions_emma.yml", emit: versions

    script:
        def base_args = (task.ext.args ?: '').toString().trim()
        def inv_flag  = (meta.invertebrates ? '--invertebrates' : '')
        def args      = [base_args, inv_flag].findAll{ it }.join(' ')
        def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
        def effective_args = [args, "--fa emma/<prefix>.fa", "--gff emma/<prefix>.gff", "--tbl emma/<prefix>.tbl", "--svg emma/<prefix>.svg", "--tempdir tempdir/", "--loglevel debug", fasta].findAll { it?.toString()?.trim() }.join(' ')

        """
            # Explicitly export the environment variable (double protection)
            export JULIA_DEPOT_PATH='/opt/julia-depot'
    
            emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
                     
            emma_prefix="${prefix}.emma"\${emma_version}""
            
            mkdir -p tempdir cds emma
            
            julia --project=/opt/Emma /opt/Emma/src/command.jl \\
                $args \\
                --fa emma/\${emma_prefix}.fa \\
                --gff emma/\${emma_prefix}.gff \\
                --tbl emma/\${emma_prefix}.tbl \\
                --svg emma/\${emma_prefix}.svg \\
                --tempdir tempdir/ \\
                --loglevel debug \\
                ${fasta} 

            cat <<-END_TOOL_PARAMS > 07_emma.tool_params_mqcrow.html
            <tr><td>EMMA</td><td><samp>${effective_args}</samp></td><td>Annotates the mitogenome assembly for ${meta.id}; adds --invertebrates when the sample metadata requires it.</td></tr>
            END_TOOL_PARAMS
            
            julia /opt/extract_proteins.jl \\
                emma/ \\
                emma/
            
            mv cds emma/
            mv proteins emma/  
            
            # Add in the sample to the file names in cds
            for file in emma/cds/*; do
                new_file="\${file%.fa}.\${emma_prefix}.fa"     
                mv "\$file" "\$new_file"
            done

            # Add in the sample to the file names in cds
            for file in emma/proteins/*; do
                new_file="\${file%.fa}.\${emma_prefix}.fa"     
                mv "\$file" "\$new_file"
            done

            cat <<-END_VERSIONS > versions_emma.yml
            "${task.process}":
                julia: \$(julia -v | sed 's/^julia version //g' )
                emma: \$(cat /opt/Emma/Project.toml | grep version)
            END_VERSIONS
    
        """

 stub:
    def prefix = task.ext.prefix ?: "${meta.mt_assembly_prefix}"
    def base_args = (task.ext.args ?: '').toString().trim()
    def inv_flag  = (meta.invertebrates ? '--invertebrates' : '')
    def args      = [base_args, inv_flag].findAll{ it }.join(' ')
    def effective_args = [args, "--fa emma/<prefix>.fa", "--gff emma/<prefix>.gff", "--tbl emma/<prefix>.tbl", "--svg emma/<prefix>.svg", "--tempdir tempdir/", "--loglevel debug", fasta].findAll { it?.toString()?.trim() }.join(' ')
    
        """
        # Create mock emma version for stub
        emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
        emma_prefix="${prefix}.emma\${emma_version}"
        
        # Create directory structure
        mkdir -p emma/cds
        
        # Create mock output files
        touch emma/\${emma_prefix}.fa
        touch emma/\${emma_prefix}.gff
        touch emma/\${emma_prefix}.tbl
        touch emma/\${emma_prefix}.svg
        
        # Create mock CDS files
        touch emma/cds/CO1_gene.\${emma_prefix}.fa
        touch emma/cds/RNR1_gene.\${emma_prefix}.fa
        touch emma/cds/RNR2_gene.\${emma_prefix}.fa

        cat <<-END_TOOL_PARAMS > 07_emma.tool_params_mqcrow.html
        <tr><td>EMMA</td><td><samp>${effective_args}</samp></td><td>Annotates the mitogenome assembly for ${meta.id}; adds --invertebrates when the sample metadata requires it.</td></tr>
        END_TOOL_PARAMS
        
        # Create mock versions file
        cat <<-END_VERSIONS > versions_emma.yml
        "${task.process}":
            julia: \$(julia -v | sed 's/^julia version //g' )
            emma: \$(cat /opt/Emma/Project.toml | grep version)
        END_VERSIONS

        """
}
