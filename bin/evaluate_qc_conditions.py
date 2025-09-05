#!/usr/bin/env python3

import csv
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Evaluate QC conditions from BLAST and annotation results')
    parser.add_argument('--blast-table', required=True, help='Path to BLAST table file')
    parser.add_argument('--annotation-csv', required=True, help='Path to annotation CSV file')
    parser.add_argument('--output-species', default='species_name.txt', help='Output file for species name')
    parser.add_argument('--output-proceed', default='proceed_qc.txt', help='Output file for proceed QC decision')
    parser.add_argument('--output-versions', default='versions.yml', help='Output file for versions')
    
    args = parser.parse_args()
    
    # Initialize variables
    proceed_qc = "false"
    species_name = "unknown"
    
    # Check blast table for Found_in_blast_YN = "Yes"
    blast_found = False
    try:
        with open(args.blast_table, 'r') as f:
            # Skip header if present
            lines = f.readlines()
            if lines:
                # Check if first line is header
                if lines[0].startswith('og_id'):
                    lines = lines[1:]
                
                for line in lines:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            found_in_blast = parts[4].strip()
                            if found_in_blast.lower() == "yes":
                                blast_found = True
                                # Extract species name from LCA_result (column 2, index 1)
                                if len(parts) >= 2:
                                    species_name = parts[1].strip()
                                break
    except Exception as e:
        print(f"Error reading blast table: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check annotation CSV for passed = "yes"
    annotation_passed = False
    try:
        with open(args.annotation_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'passed' in row and row['passed'].lower() == "yes":
                    annotation_passed = True
                    break
    except Exception as e:
        print(f"Error reading annotation CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Both conditions must be satisfied
    if blast_found and annotation_passed:
        proceed_qc = "true"
    
    # Write results to files for Nextflow
    with open(args.output_species, 'w') as f:
        f.write(species_name)
    
    with open(args.output_proceed, 'w') as f:
        f.write(proceed_qc)
    
    print(f"Blast condition met: {blast_found}")
    print(f"Annotation condition met: {annotation_passed}")
    print(f"Species name: {species_name}")
    print(f"Proceed with QC: {proceed_qc}")
    
    # Create versions file
    with open(args.output_versions, 'w') as f:
        f.write(f'"{os.environ.get("NF_TASK_PROCESS", "EVALUATE_QC_CONDITIONS")}":\n')
        f.write(f'    python: "{sys.version.split()[0]}"\n')

if __name__ == "__main__":
    main()