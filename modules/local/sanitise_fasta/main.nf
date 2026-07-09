process SANITISE_FASTA {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('output/*.{fa,fasta,fna}'), emit: fasta

    when:
    !params.skip_mitogenome_annotation

    script:
    // Derive the prefix from the input fasta basename, NOT meta.mt_assembly_prefix:
    // the reseed pass carries its "reseed" suffix only on the filename (meta keeps the
    // first-pass "getorg<ver>"), so keying off meta mislabels a reseed concat as
    // "getorg<ver>_concat" instead of "getorg<ver>reseed_concat". Strip only the FASTA
    // extension (simpleName would collapse every dot-field down to the OG id). Fall
    // back to mt_assembly_prefix / id only if the basename somehow comes back empty.
    def input_name = fasta.name
    def base = input_name.replaceAll(/\.(fa|fasta|fna)$/, '')
    def prefix = base ?: (meta.containsKey('mt_assembly_prefix') && meta.mt_assembly_prefix ? meta.mt_assembly_prefix : meta.id)
    def concat_name = prefix + '_concat'
    def extension = input_name.tokenize('.').last()
    def output_concat = "output/${concat_name}.${extension}"
    def output_single = "output/${input_name}"
    
    """
    set -euo pipefail
    
    # Create output directory
    mkdir -p output
    
    # Count contigs
    n=\$(grep -c '^>' ${fasta} || true)
    
    if [ "\$n" -gt 1 ]; then
        # Multiple contigs: concatenate all sequences into a single record
        echo "Multiple contigs detected (\$n), concatenating..."
        awk 'BEGIN{ORS=""} /^>/ {next} {gsub(/[ \\t\\r\\n]/,""); print \$0} END{print "\\n"}' ${fasta} > seq.tmp
        printf ">%s\\n" "${concat_name}" > ${output_concat}
        cat seq.tmp >> ${output_concat}
        rm -f seq.tmp
        echo "Created concatenated file: ${output_concat}"
    else
        # Single contig: copy to output directory with original name
        echo "Single contig detected, copying to output directory..."
        cp ${fasta} ${output_single}
        echo "Created output file: ${output_single}"
    fi
    
    echo "Output directory contents:"
    ls -la output/
    """
}