process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/889a182b8066804f4799f3808a5813ad601381a8a0e3baa4ab8d73e739b97001/data' :
        'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690' }"

    input:
    tuple val(meta), path(reads)
    path  adapter_fasta

    output:
    tuple val(meta), path('*.fastp.fastq.gz'), emit: reads
    tuple val(meta), path('*.fastp.json')    , emit: json
    tuple val(meta), path('*.fastp.html')    , emit: html
    tuple val(meta), path('*.fastp.log')     , emit: log
    tuple val(meta), path('04_fastp.tool_params_mqcrow.html'), emit: tool_params
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta.sequencing_type}.${meta.date}"
    def adapter_arg = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ''
    def effective_args = [args, adapter_arg, "--thread ${task.cpus}"].findAll { it?.toString()?.trim() }.join(' ')
    def note        = 'Adapter + quality trimming of raw HiC reads (poly-G removal) before assembly.'
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}.fastp.fastq.gz \\
            ${adapter_arg} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --report_title="${prefix} fastp" \\
            --thread ${task.cpus} \\
            ${args} \\
            2>&1 | tee ${prefix}.fastp.log

        cat <<-END_TOOL_PARAMS > 04_fastp.tool_params_mqcrow.html
        <tr><td>FASTP</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}.R1.fastp.fastq.gz \\
            --in2 ${reads[1]} \\
            --out2 ${prefix}.R2.fastp.fastq.gz \\
            ${adapter_arg} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --report_title="${prefix} fastp" \\
            --thread ${task.cpus} \\
            ${args} \\
            2>&1 | tee ${prefix}.fastp.log

        cat <<-END_TOOL_PARAMS > 04_fastp.tool_params_mqcrow.html
        <tr><td>FASTP</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
        END_TOOL_PARAMS

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta.sequencing_type}.${meta.date}"
    def adapter_arg = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ''
    def effective_args = [args, adapter_arg, "--thread ${task.cpus}"].findAll { it?.toString()?.trim() }.join(' ')
    def note        = 'Adapter + quality trimming of raw HiC reads (poly-G removal) before assembly.'
    def touch_reads = meta.single_end ?
        "echo '' | gzip > ${prefix}.fastp.fastq.gz" :
        "echo '' | gzip > ${prefix}.R1.fastp.fastq.gz ; echo '' | gzip > ${prefix}.R2.fastp.fastq.gz"
    """
    ${touch_reads}
    touch ${prefix}.fastp.json
    touch ${prefix}.fastp.html
    touch ${prefix}.fastp.log

    cat <<-END_TOOL_PARAMS > 04_fastp.tool_params_mqcrow.html
    <tr><td>FASTP</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
