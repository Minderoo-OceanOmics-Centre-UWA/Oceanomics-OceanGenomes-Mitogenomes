#!/usr/bin/env python3

import csv
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Evaluate QC conditions from BLAST and annotation results')
    parser.add_argument('--blast-table', required=True, help='Path to BLAST table file')
    parser.add_argument('--annotation-csv', required=True, help='Path to annotation CSV file')
    parser.add_argument('--circularity-check', help='Path to circularity-check TSV (length/repeat anomaly gate)')
    parser.add_argument('--circular', default='null',
                        help='meta.circular from the assembler (GetOrganelle log verdict): '
                             'true / false / null. MitoHiFi leaves this null and the verdict '
                             'is taken from --circularity-check (final_verdict_circular) instead.')
    parser.add_argument('--output-species', default='species_name.txt', help='Output file for species name')
    parser.add_argument('--output-proceed', default='proceed_qc.txt', help='Output file for proceed QC decision')
    parser.add_argument('--output-circular', default='circular.txt',
                        help='Output file for the resolved circular verdict (true/false). Consumed '
                             'downstream to set the table2asn topology/completeness modifiers.')
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
                                # Extract species name from nom_species_id (column 3, index 1)
                                if len(parts) >= 3:
                                    species_name = parts[2].strip()
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
    
    # Check the circularity-check sidecar for a length / tandem-repeat anomaly. A
    # concatemer (tandem genome duplication), control_region_repeat (D-loop VNTR)
    # or unresolved over-length assembly is not submission-ready and must not
    # progress to QC until curated. Absent/empty file (e.g. GetOrganelle samples or
    # the placeholder) -> no anomaly, no block.
    assembly_anomaly = False
    anomaly_type = "none"
    # Circularity verdict from the (MitoHiFi) sidecar, read alongside the anomaly
    # check so the file is only opened once. None = no information in this file.
    circ_from_file = None
    if args.circularity_check and os.path.exists(args.circularity_check):
        try:
            with open(args.circularity_check, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    at = (row.get('anomaly_type') or '').strip().lower()
                    la = (row.get('length_anomaly') or '').strip().lower()
                    if at not in ('', 'none', 'na'):
                        assembly_anomaly = True
                        anomaly_type = at
                    elif la == 'yes':
                        assembly_anomaly = True
                        anomaly_type = 'length_anomaly'
                    fv = (row.get('final_verdict_circular') or '').strip().lower()
                    if fv == 'true':
                        circ_from_file = True
                    elif fv == 'false':
                        circ_from_file = False
                    # 'na'/'' (GetOrganelle placeholder) -> leave as None (no info)
                    break
        except Exception as e:
            print(f"Error reading circularity-check: {e}", file=sys.stderr)

    # Resolve the circular verdict. GetOrganelle reports it via meta.circular
    # (--circular); MitoHiFi via the sidecar's final_verdict_circular. Exactly one
    # source is informative per sample; the other is null/NA. "unknown" (e.g.
    # precomputed reruns, or a log parse that failed) is treated conservatively:
    # it is NOT a circular genome, so the topology defaults to linear, but it is
    # also NOT a confirmed non-circular assembly, so it does not by itself block QC.
    meta_circ = (args.circular or 'null').strip().lower()
    circular_true = (meta_circ == 'true') or (circ_from_file is True)
    circular_false = (meta_circ == 'false') or (circ_from_file is False)
    circular_known_false = circular_false and not circular_true
    circular_out = "true" if circular_true else "false"

    # All conditions must be satisfied: species found, annotation passed, no
    # length/repeat anomaly, and the assembly is not a confirmed non-circular
    # molecule. A non-circular mitogenome is not submission-ready, so it must not
    # pass QC (it would otherwise be exported with a false circular topology).
    if blast_found and annotation_passed and not assembly_anomaly and not circular_known_false:
        proceed_qc = "true"

    # Write results to files for Nextflow
    with open(args.output_species, 'w') as f:
        f.write(species_name)

    with open(args.output_proceed, 'w') as f:
        f.write(proceed_qc)

    with open(args.output_circular, 'w') as f:
        f.write(circular_out)

    print(f"Blast condition met: {blast_found}")
    print(f"Annotation condition met: {annotation_passed}")
    print(f"Assembly anomaly: {assembly_anomaly} ({anomaly_type})")
    print(f"Circular verdict: {circular_out} (meta={meta_circ}, file={circ_from_file}, known_non_circular={circular_known_false})")
    print(f"Species name: {species_name}")
    print(f"Proceed with QC: {proceed_qc}")
    
    # Create versions file
    with open(args.output_versions, 'w') as f:
        f.write(f'"{os.environ.get("NF_TASK_PROCESS", "EVALUATE_QC_CONDITIONS")}":\n')
        f.write(f'    python: "{sys.version.split()[0]}"\n')

if __name__ == "__main__":
    main()