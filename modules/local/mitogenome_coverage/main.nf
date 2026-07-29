// Uniform, cross-platform mitogenome read depth.
//
// Replaces three mutually incomparable "coverage" numbers with one definition:
// mean per-base depth of the sample's own reads remapped to the assembly that
// actually goes to annotation. Previously GetOrganelle reported k-mer coverage
// off its assembly graph (~0.2x true depth, and measured on the reduced read set
// its --reduce-reads-for-coverage default selects), MitoHiFi reported depth of
// only the reads it had recruited by mapping to a related-species reference (so a
// divergent reference silently depressed the number), and Oatk reported nothing
// at all. None of those could be compared to each other.
//
// Remapping to the sample's OWN assembly rather than to a reference is the point:
// a genuinely divergent mitogenome is measured just as accurately as a
// well-referenced one, which is exactly the bias that made the MitoHiFi number
// untrustworthy.
//
// Runs on the SANITISE_FASTA output, i.e. the exact molecule that reaches
// annotation and GenBank. That placement is deliberate: it is post-collapse (a
// concatemer would otherwise halve the reported depth), post-reseed (only the
// variant that won is worth measuring), and post-concatenation (SANITISE_FASTA
// rewrites a multi-contig assembly into a single _concat record, so the
// pre-sanitise FASTA is not the molecule anyone will look at). It also means no
// remap is ever spent on an assembly that failed or fell below the length floor.
//
// Reuses the pinned MitoHiFi image, which ships the minimap2 + python3 this
// needs, exactly as REFERENCE_RANK already does for short-read candidate ranking.
process MITOGENOME_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    container "${params.mitohifi_container}"

    input:
    tuple val(meta), path(fasta), path(reads)

    output:
    tuple val(meta), path("${depth_prefix}.mito_depth.tsv"), emit: depth
    tuple val(meta), path("06b_mitogenome_coverage.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // HiFi needs map-hifi; Illumina/HiC are short reads and need sr. Getting this
    // wrong would map essentially nothing and report a depth near zero.
    def preset = meta.sequencing_type == 'hifi' ? 'map-hifi' : 'sr'
    def minId = meta.sequencing_type == 'hifi'
        ? (params.mitogenome_depth_min_identity_hifi ?: 0.99)
        : (params.mitogenome_depth_min_identity_sr ?: 0.95)
    def minFrac = meta.sequencing_type == 'hifi'
        ? (params.mitogenome_depth_min_aligned_frac_hifi ?: 0.70)
        : (params.mitogenome_depth_min_aligned_frac_sr ?: 0.80)
    def fraction = params.mitogenome_depth_subsample_fraction ?: 0
    def circular = meta.circular == null ? 'null' : meta.circular.toString()
    // Name the output after the FASTA actually measured, NOT meta.mt_assembly_prefix:
    // the reseed / _rgj suffix lives only on the filename (meta keeps the first-pass
    // prefix), and the assembly summary groups runs by filename. Strip a trailing
    // _collapsed / _concat because neither has a summary run of its own -- leaving
    // either on would manufacture a phantom row carrying a depth but no assembly.
    depth_prefix = fasta.baseName.replaceAll(/_(collapsed|concat)$/, '')
    def effective_args = [
        "minimap2 -ax ${preset} --secondary=no",
        "--min-identity ${minId}",
        "--min-aligned-frac ${minFrac}",
        "--circular ${circular}",
        fraction ? "--subsample-fraction ${fraction}" : '',
        args
    ].findAll { it?.toString()?.trim() }.join(' ')
    """
    mito_depth.py \\
        --fasta ${fasta} \\
        --reads ${reads} \\
        --sample ${meta.mt_assembly_prefix ?: meta.id} \\
        --out ${depth_prefix}.mito_depth.tsv \\
        --preset ${preset} \\
        --sequencing-type ${meta.sequencing_type ?: ''} \\
        --circular ${circular} \\
        --min-identity ${minId} \\
        --min-aligned-frac ${minFrac} \\
        --subsample-fraction ${fraction} \\
        --threads ${task.cpus} \\
        --workdir . \\
        ${args}

    cat <<-END_TOOL_PARAMS > 06b_mitogenome_coverage.tool_params_mqcrow.html
    <tr><td>Mitogenome Coverage</td><td><samp>${effective_args}</samp></td><td>Remaps the full read set for ${meta.id} to its final assembly and reports mean per-base depth, folded on a doubled reference when the molecule is circular.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>/dev/null)
        python: \$(python3 --version | sed 's/^Python //')
    END_VERSIONS
    """

    stub:
    def preset = meta.sequencing_type == 'hifi' ? 'map-hifi' : 'sr'
    depth_prefix = fasta.baseName.replaceAll(/_(collapsed|concat)$/, '')
    """
    printf 'sample\\ttarget_fasta\\ttarget_length_bp\\tn_contigs\\tcircular_doubled\\tpreset\\tsequencing_type\\tmean_depth\\tmedian_depth\\tsd_depth\\tdepth_cv\\tp10_depth\\tmin_depth\\tbreadth_1x\\tbreadth_10x\\tbreadth_20x\\tmito_mapped_reads\\treads_fail_identity\\treads_fail_clip\\tsupplementary_dropped\\ttotal_reads\\tmito_read_fraction\\tmean_identity\\tmin_identity\\tmin_aligned_frac\\tsubsampled\\tsubsample_fraction\\tscale_factor\\tdepth_method\\tmean_coverage\\tcoverage_cv\\n%s\\t%s\\t16500\\t1\\tfalse\\t%s\\t%s\\t100\\t100\\t10\\t0.1\\t90\\t80\\t1\\t1\\t1\\t1000\\t0\\t0\\t0\\t10000\\t0.1\\t0.999\\t0.95\\t0.8\\tfalse\\tNA\\t1\\tremap_full_v1\\t100\\t0.1\\n' \\
        '${meta.mt_assembly_prefix ?: meta.id}' '${fasta}' '${preset}' '${meta.sequencing_type ?: ''}' \\
        > ${depth_prefix}.mito_depth.tsv

    cat <<-END_TOOL_PARAMS > 06b_mitogenome_coverage.tool_params_mqcrow.html
    <tr><td>Mitogenome Coverage</td><td><samp>mito_depth.py (stub)</samp></td><td>Remaps the full read set for ${meta.id} to its final assembly.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "stub"
        python: "stub"
    END_VERSIONS
    """
}
