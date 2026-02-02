#!/usr/bin/env python3
# singularity run $SING/psycopg2:0.1.sif python3 create_samplesheet.py --output samplesheet.csv --sql-config "/home/tpeirce/postgresql_details/oceanomics.cfg" --input-files /scratch/pawsey0964/tpeirce/BIGRUN/mitogenomes/NOVA_251215_AD/fastp/*/fastp/*.fastq.gz
import argparse
import csv
import glob
import os
import re
import sys
import configparser
from datetime import datetime, timedelta
from pathlib import Path

import psycopg2


def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"Config file '{config_file}' does not exist.")

    config = configparser.ConfigParser()
    config.read(config_file)

    if not config.has_section('postgres'):
        raise ValueError("Missing [postgres] section in config file.")

    required_keys = ['dbname', 'user', 'password', 'host', 'port']
    for key in required_keys:
        if not config.has_option('postgres', key):
            raise ValueError(f"Missing '{key}' in [postgres] section of config file.")

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')
    }


def extract_sample_info(filename):
    """
    Extract sample information from filename based on three naming patterns:
    1. OG820M-1_HICL_S4_L001_R1_001.fastq.gz -> OG820M-1_HICL (before second _)
    2. OG785_m84154_241004_105305_s3.hifi_reads.bc2068.filt.fastq.gz -> OG785 (before first _)
    3. OG764.ilmn.240716.R1.fastq.gz -> OG764 (before first .)
    """
    basename = Path(filename).name

    if '_' in basename:
        parts = basename.split('_')
        if len(parts) >= 2 and re.match(r'^HIC[A-Za-z0-9]*$', parts[1]):
            return f"{parts[0]}_{parts[1]}"
        return parts[0]
    elif '.' in basename:
        return basename.split('.')[0]
    return basename.replace('.fastq.gz', '').replace('.fq.gz', '').replace('.fastq', '').replace('.fq', '')


def determine_sequencing_type(filename):
    basename = Path(filename).name.lower()

    if 'hifi_reads' in basename or 'hifi.reads' in basename:
        return 'hifi'
    elif 'hicl' in basename:
        return 'hic'
    elif 'ilmn' in basename:
        return 'ilmn'
    return 'unknown'


def determine_file_type(filename):
    basename = Path(filename).name.lower()

    if any(pattern in basename for pattern in ['_r1_', '_r1.', '.r1.', '_1_', '_1.']):
        return 'R1'
    elif any(pattern in basename for pattern in ['_r2_', '_r2.', '.r2.', '_2_', '_2.']):
        return 'R2'
    elif any(pattern in basename for pattern in ['hifi_reads', 'hifi.reads', '.hifi']):
        return 'single'
    return 'single'


def group_files_by_sample(file_list):
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


def pick_primary_seq_type(seq_types, sample_id):
    if len(seq_types) == 1:
        return list(seq_types)[0]

    if 'hifi' in seq_types:
        primary_seq_type = 'hifi'
    elif 'hic' in seq_types:
        primary_seq_type = 'hic'
    elif 'ilmn' in seq_types:
        primary_seq_type = 'ilmn'
    else:
        primary_seq_type = 'mixed'

    print(
        f"Warning: Sample {sample_id} has multiple sequencing types: {seq_types}. "
        f"Using: {primary_seq_type}",
        file=sys.stderr
    )
    return primary_seq_type


def clean_sample_id(sample_id):
    return re.sub(r'^(OG\\d+).*', r'\\1', sample_id)


def extract_hifi_completion_date(filename):
    parts = Path(filename).name.split('_')
    return parts[2] if len(parts) >= 3 else 'unknown'


def extract_ilmn_date(filename):
    parts = Path(filename).name.split('.')
    return parts[2] if len(parts) >= 3 else 'unknown'


def extract_ilmn_prefix(filename, sample_id, sequencing_type, date):
    parts = Path(filename).name.split('.')
    if len(parts) >= 3:
        return f"{parts[0]}.{parts[1]}.{parts[2]}"
    return f"{sample_id}.{sequencing_type}.{date}"


def query_hifi_date(cursor, sample_id, completion_date_str):
    if cursor is None:
        return "unknown"
    if completion_date_str and completion_date_str != "unknown":
        try:
            completion_date = datetime.strptime(completion_date_str, '%y%m%d').date()
            search_start = completion_date - timedelta(days=30)
            cursor.execute(
                """
                SELECT seq_date
                FROM sequencing
                WHERE og_id = %s
                  AND seq_date >= %s
                  AND seq_date <= %s
                ORDER BY seq_date DESC
                LIMIT 1
                """,
                (sample_id, search_start.strftime('%y%m%d'), completion_date.strftime('%y%m%d'))
            )
        except ValueError as exc:
            print(f"Error parsing completion date '{completion_date_str}': {exc}", file=sys.stderr)
            cursor.execute(
                "SELECT seq_date FROM sequencing WHERE og_id = %s ORDER BY seq_date DESC LIMIT 1",
                (sample_id,)
            )
    else:
        return "unknown"

    result = cursor.fetchone()
    return str(result[0]) if result else "unknown"


def query_hic_date(cursor, sample_id):
    if cursor is None:
        return "unknown"
    cursor.execute(
        "SELECT seq_date FROM sequencing WHERE hic_library_tube_id = %s",
        (sample_id,)
    )
    result = cursor.fetchone()
    if not result:
        return "unknown"
    value = result[0]
    return value.strftime('%y%m%d') if hasattr(value, 'strftime') else str(value)


def query_nominal_species_id(cursor, sample_id):
    if cursor is None:
        return "unknown"
    cursor.execute(
        "SELECT nominal_species_id FROM sample WHERE og_id = %s",
        (sample_id,)
    )
    result = cursor.fetchone()
    if not result:
        return "unknown"
    value = result[0]
    return value.strftime('%y%m%d') if hasattr(value, 'strftime') else str(value)


def parse_args():
    parser = argparse.ArgumentParser(description="Create enriched samplesheet from staged FASTQ files.")
    parser.add_argument("--output", required=True, help="Output samplesheet CSV path.")
    parser.add_argument("--sql-config", required=False, default=None, help="Optional path to SQL config file.")
    parser.add_argument("--input-files", nargs="*", default=None, help="Optional explicit list of input files.")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.input_files:
        file_list = [os.path.abspath(f) for f in args.input_files]
    else:
        file_list = []
        for ext in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
            file_list.extend(glob.glob(ext))
        file_list = [os.path.abspath(f) for f in file_list]

    header = [
        'sample',
        'sequencing_type',
        'single_end',
        'original_id',
        'completion_date',
        'date',
        'assembly_prefix',
        'nominal_species_id',
        'fastq_1',
        'fastq_2'
    ]

    if not file_list:
        print("Warning: No FASTQ files found in staged directory", file=sys.stderr)
        with open(args.output, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
        return 0

    samples = group_files_by_sample(file_list)

    cursor = None
    conn = None
    if args.sql_config:
        try:
            db_params = load_db_config(args.sql_config)
            conn = psycopg2.connect(**db_params)
            cursor = conn.cursor()
        except Exception as exc:
            print(f"Error connecting to database: {exc}", file=sys.stderr)

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for sample_id, files in sorted(samples.items()):
            seq_types = files['sequencing_types']
            sequencing_type = pick_primary_seq_type(seq_types, sample_id)

            original_id = sample_id
            cleaned_id = clean_sample_id(sample_id) if sequencing_type in {'hifi', 'hic'} else sample_id
            completion_date = ''

            first_file = None
            if files['single']:
                first_file = files['single'][0]
            elif files['R1']:
                first_file = files['R1'][0]
            elif files['R2']:
                first_file = files['R2'][0]

            if sequencing_type == 'hifi':
                completion_date = extract_hifi_completion_date(first_file) if first_file else 'unknown'
                date = query_hifi_date(cursor, cleaned_id, completion_date)
            elif sequencing_type == 'hic':
                date = query_hic_date(cursor, original_id)
            elif sequencing_type == 'ilmn':
                date = extract_ilmn_date(first_file) if first_file else 'unknown'
            else:
                date = 'unknown'

            if sequencing_type == 'ilmn':
                assembly_prefix = extract_ilmn_prefix(first_file, cleaned_id, sequencing_type, date) if first_file else ''
            else:
                assembly_prefix = f"{cleaned_id}.{sequencing_type}.{date}"

            nominal_species_id = query_nominal_species_id(cursor, cleaned_id)

            def write_row(r1, r2, single_end):
                writer.writerow([
                    cleaned_id,
                    sequencing_type,
                    str(single_end).lower(),
                    original_id,
                    completion_date,
                    date,
                    assembly_prefix,
                    nominal_species_id,
                    r1,
                    r2
                ])

            if files['R1'] and files['R2']:
                r1_files = sorted(files['R1'])
                r2_files = sorted(files['R2'])
                max_files = max(len(r1_files), len(r2_files))
                for i in range(max_files):
                    r1_file = r1_files[i] if i < len(r1_files) else ''
                    r2_file = r2_files[i] if i < len(r2_files) else ''
                    if r1_file or r2_file:
                        write_row(r1_file, r2_file, False)
            elif files['single']:
                for single_file in sorted(files['single']):
                    write_row(single_file, '', True)
            elif files['R1']:
                for r1_file in sorted(files['R1']):
                    write_row(r1_file, '', True)
            elif files['R2']:
                for r2_file in sorted(files['R2']):
                    write_row(r2_file, '', True)

    if cursor is not None:
        cursor.close()
    if conn is not None:
        conn.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
