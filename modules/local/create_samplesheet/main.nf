process CREATE_SAMPLESHEET {
    tag "Creating samplesheet from ${input_path}"
    
    conda "conda-forge::python=3.9"
    
    input:
    path(input_files)
    val output_name
    
    output:
    path "${output_name}", emit: samplesheet
    path "versions.yml", emit: versions
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    import glob
    import re
    import csv
    from pathlib import Path
    
    def extract_sample_info(filename):
        \"\"\"
        Extract sample information from filename based on three naming patterns:
        1. OG820M-1_HICL_S4_L001_R1_001.fastq.gz -> OG820M-1_HICL (before second _)
        2. OG785_m84154_241004_105305_s3.hifi_reads.bc2068.filt.fastq.gz -> OG785 (before first _)
        3. OG764.ilmn.240716.R1.fastq.gz -> OG764 (before first .)
        \"\"\"
        basename = Path(filename).name
        
        # Underscore-based names (patterns 1 & 2)
        if '_' in basename:
            parts = basename.split('_')
            # If second token looks like a Hi-C tag, return first two parts
            if len(parts) >= 2 and re.match(r'^HIC[A-Za-z0-9]*\$', parts[1]):
                return f"{parts[0]}_{parts[1]}"
            # Otherwise just the first token
            return parts[0]
        
        # If no underscore, check for dot (format 3)
        elif '.' in basename:
            # Extract everything before the first dot
            sample_id = basename.split('.')[0]
            return sample_id
        
        else:
            # Fallback: return the whole basename without extension
            return basename.replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')
    
    def determine_sequencing_type(filename):
        \"\"\"
        Determine sequencing type based on filename patterns:
        - HiC: contains 'HICL'
        - HiFi: contains 'hifi_reads' or 'hifi.reads'
        - Illumina: contains 'ilmn'
        \"\"\"
        basename = Path(filename).name.lower()
        
        if 'hicl' in basename:
            return 'hic'
        elif 'hifi_reads' in basename or 'hifi.reads' in basename:
            return 'hifi'
        elif 'ilmn' in basename:
            return 'ilmn'
        else:
            # Try to infer from other patterns if the main indicators aren't present
            if any(pattern in basename for pattern in ['_r1_', '_r2_', '.r1.', '.r2.']):
                return 'ilmn'  # Assume paired-end is Illumina
            else:
                return 'unknown'
    
    def determine_file_type(filename):
        \"\"\"Determine if file is R1, R2, or single-end based on filename patterns\"\"\"
        basename = Path(filename).name.lower()
        
        # Check for R1 patterns
        if any(pattern in basename for pattern in ['_r1_', '_r1.', '.r1.', '_1_', '_1.']):
            return 'R1'
        # Check for R2 patterns  
        elif any(pattern in basename for pattern in ['_r2_', '_r2.', '.r2.', '_2_', '_2.']):
            return 'R2'
        # Check for HiFi reads or other single-end indicators
        elif any(pattern in basename for pattern in ['hifi_reads', 'hifi.reads', '.hifi', 'hicl']):
            return 'single'
        else:
            return 'single'
    
    def group_files_by_sample(file_list):
        \"\"\"Group files by sample ID and organize by type and sequencing platform\"\"\"
        samples = {}
        
        for filepath in file_list:
            sample_id = extract_sample_info(filepath)
            file_type = determine_file_type(filepath)
            seq_type = determine_sequencing_type(filepath)
            
            if sample_id not in samples:
                samples[sample_id] = {
                    'R1': [], 
                    'R2': [], 
                    'single': [],
                    'sequencing_types': set()
                }
            
            samples[sample_id][file_type].append(filepath)
            samples[sample_id]['sequencing_types'].add(seq_type)
        
        return samples
    
    def create_samplesheet(samples, output_file):
        \"\"\"Create samplesheet with proper format for mixed file types including sequencing type\"\"\"
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header with sequencing_type column
            writer.writerow(['sample', 'fastq_1', 'fastq_2', 'sequencing_type'])
            
            for sample_id, files in sorted(samples.items()):
                # Determine the primary sequencing type for this sample
                seq_types = files['sequencing_types']
                if len(seq_types) == 1:
                    primary_seq_type = list(seq_types)[0]
                else:
                    # If multiple types, prioritize in order: hifi, hic, ilmn, unknown
                    if 'hifi' in seq_types:
                        primary_seq_type = 'hifi'
                    elif 'hic' in seq_types:
                        primary_seq_type = 'hic'
                    elif 'ilmn' in seq_types:
                        primary_seq_type = 'ilmn'
                    else:
                        primary_seq_type = 'mixed'
                    
                    print(f"Warning: Sample {sample_id} has multiple sequencing types: {seq_types}. Using: {primary_seq_type}")
                
                # Handle paired-end files (R1 and R2)
                if files['R1'] and files['R2']:
                    # Sort files to ensure consistent pairing
                    r1_files = sorted(files['R1'])
                    r2_files = sorted(files['R2'])
                    
                    # Pair R1 and R2 files - if multiple files per sample, they will auto-concatenate
                    max_files = max(len(r1_files), len(r2_files))
                    for i in range(max_files):
                        r1_file = r1_files[i] if i < len(r1_files) else ''
                        r2_file = r2_files[i] if i < len(r2_files) else ''
                        
                        if r1_file or r2_file:
                            writer.writerow([sample_id, r1_file, r2_file, primary_seq_type])
                    
                    # Warn if R1 and R2 counts don't match
                    if len(r1_files) != len(r2_files):
                        print(f"Warning: Unequal R1 ({len(r1_files)}) and R2 ({len(r2_files)}) files for sample {sample_id}")
                
                # Handle single-end files (including HiFi reads and HiC)
                elif files['single']:
                    for single_file in sorted(files['single']):
                        writer.writerow([sample_id, single_file, '', primary_seq_type])
                
                # Handle cases with only R1 or only R2 (treat as single-end)
                elif files['R1']:
                    for r1_file in sorted(files['R1']):
                        writer.writerow([sample_id, r1_file, '', primary_seq_type])
                elif files['R2']:
                    for r2_file in sorted(files['R2']):
                        writer.writerow([sample_id, r2_file, '', primary_seq_type])
    
    # Main execution - Get files from current directory (Nextflow stages them here)
    output_filename = "${output_name}"
    
    # Get all FASTQ files in the current directory (Nextflow stages input files here)
    file_list = []
    for ext in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
        file_list.extend(glob.glob(ext))
    
    # Convert to absolute paths
    file_list = [os.path.abspath(f) for f in file_list]
    
    if not file_list:
        print("Warning: No FASTQ files found in staged directory")
        # Create empty samplesheet
        with open(output_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['sample', 'fastq_1', 'fastq_2', 'sequencing_type'])
    else:
        print(f"Found {len(file_list)} FASTQ files")
        for f in file_list:
            print(f"  - {f}")
        
        # Group files by sample
        samples = group_files_by_sample(file_list)
        
        print(f"Identified {len(samples)} unique samples:")
        for sample_id, files in samples.items():
            r1_count = len(files['R1'])
            r2_count = len(files['R2'])
            single_count = len(files['single'])
            seq_types = ', '.join(sorted(files['sequencing_types']))
            print(f"  {sample_id}: R1={r1_count}, R2={r2_count}, single={single_count}, types={seq_types}")
        
        # Create the samplesheet
        create_samplesheet(samples, output_filename)
        
        print(f"Samplesheet created: {output_filename}")
    
    # Create versions file
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "3.9"\\n')
    """
}