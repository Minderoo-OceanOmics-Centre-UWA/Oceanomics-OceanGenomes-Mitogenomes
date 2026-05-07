process MULTIQC_PER_SAMPLE {
    label 'process_medium'

    conda "${projectDir}/modules/nf-core/multiqc/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_1' :
        'biocontainers/multiqc:1.30--pyhdfd78af_1' }"

    input:
    path  multiqc_files, stageAs: "multiqc_per_sample_inputs/?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    path "sample_multiqc_html/*/*/*_multiqc_report.html", emit: reports
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    def multiqc_arg_chunks = groovy.json.JsonOutput.toJson([args, config, extra_config, logo, replace, samples])
    """
    cat > multiqc_per_sample_args.json <<'EOF'
${multiqc_arg_chunks}
EOF

    python3 "$projectDir/bin/multiqc_per_sample.py" \\
        --input-root multiqc_per_sample_inputs \\
        --sample-input-root per_sample_inputs \\
        --output-root per_sample \\
        --html-output-root sample_multiqc_html \\
        --multiqc-arg-chunks-file multiqc_per_sample_args.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sample_multiqc_html/STUB/STUB
    touch sample_multiqc_html/STUB/STUB/STUB_multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "1.30"
    END_VERSIONS
    """
}
