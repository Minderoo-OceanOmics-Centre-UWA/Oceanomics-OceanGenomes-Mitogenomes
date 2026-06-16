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
    // Per-gene CDS outputs are optional: a partial/divergent mitogenome may not
    // yield every gene. Without optional, a missing CO1/RNR1/RNR2 turns a
    // successful annotation into a failed task and the sample is dropped.
    tuple val(meta), path("annotation/cds/*CO1*.fa"),  emit: co1_sequences, optional: true
    tuple val(meta), path("annotation/cds/*RNR1*.fa"),  emit: s12_sequences, optional: true
    tuple val(meta), path("annotation/cds/*RNR2*.fa"),  emit: s16_sequences, optional: true
    tuple val(meta), path("annotation/*"), emit: results
    tuple val(meta), path("07_emma.tool_params_mqcrow.html"), emit: tool_params
    path "versions_emma.yml", emit: versions

    script:
        def base_args = (task.ext.args ?: '').toString().trim()
        def inv_flag  = (meta.invertebrates ? '--invertebrates' : '')
        def args      = [base_args, inv_flag].findAll{ it }.join(' ')
        def prefix = task.ext.prefix ?: meta.mt_assembly_prefix
        def effective_args = [args, "--fa annotation/<prefix>.fa", "--gff annotation/<prefix>.gff", "--tbl annotation/<prefix>.tbl", "--svg annotation/<prefix>.svg", "--tempdir tempdir/", "--loglevel debug", fasta].findAll { it?.toString()?.trim() }.join(' ')

        """
            # Explicitly export the environment variable (double protection)
            export JULIA_DEPOT_PATH='/opt/julia-depot'
    
            emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
                     
            emma_prefix="${prefix}.emma"\${emma_version}""
            
            mkdir -p tempdir cds annotation
            
            julia --project=/opt/Emma /opt/Emma/src/command.jl \\
                $args \\
                --fa annotation/\${emma_prefix}.fa \\
                --gff annotation/\${emma_prefix}.gff \\
                --tbl annotation/\${emma_prefix}.tbl \\
                --svg annotation/\${emma_prefix}.svg \\
                --tempdir tempdir/ \\
                --loglevel debug \\
                ${fasta} 

            cat <<-END_TOOL_PARAMS > 07_emma.tool_params_mqcrow.html
            <tr><td>EMMA</td><td><samp>${effective_args}</samp></td><td>Annotates the mitogenome assembly for ${meta.id}; adds --invertebrates when the sample metadata requires it.</td></tr>
            END_TOOL_PARAMS
            
            julia /opt/extract_proteins.jl \\
                annotation/ \\
                annotation/
            
            mv cds annotation/
            mv proteins annotation/  
            
            # Add in the sample to the file names in cds
            for file in annotation/cds/*; do
                new_file="\${file%.fa}.\${emma_prefix}.fa"     
                mv "\$file" "\$new_file"
            done

            # Add in the sample to the file names in cds
            for file in annotation/proteins/*; do
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
    def effective_args = [args, "--fa annotation/<prefix>.fa", "--gff annotation/<prefix>.gff", "--tbl annotation/<prefix>.tbl", "--svg annotation/<prefix>.svg", "--tempdir tempdir/", "--loglevel debug", fasta].findAll { it?.toString()?.trim() }.join(' ')
    
        """
        # Create mock emma version for stub
        emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
        emma_prefix="${prefix}.emma\${emma_version}"
        
        # Create directory structure
        mkdir -p annotation/cds
        
        # Create mock output files
        touch annotation/\${emma_prefix}.fa
        touch annotation/\${emma_prefix}.gff
        touch annotation/\${emma_prefix}.tbl
        touch annotation/\${emma_prefix}.svg
        
        # Create mock CDS files
        touch annotation/cds/CO1_gene.\${emma_prefix}.fa
        touch annotation/cds/RNR1_gene.\${emma_prefix}.fa
        touch annotation/cds/RNR2_gene.\${emma_prefix}.fa

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
