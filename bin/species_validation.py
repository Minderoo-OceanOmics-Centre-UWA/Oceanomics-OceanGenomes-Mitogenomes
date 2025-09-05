#!/usr/bin/env python3
import psycopg2
import csv
import configparser
import sys

def load_db_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')    
    }

def concatenate_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                outfile.writelines(infile)

def normalise_name(name):
    return name.strip().lower().replace('_', ' ')

def get_species_for_ogid(db_params, og_id):
    query = "SELECT nominal_species_id FROM sample WHERE og_id = %s"
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            cur.execute(query, (og_id,))
            result = cur.fetchone()
            return result[0] if result else None
    except Exception as e:
        print(f"[ERROR] Database access failed: {e}")
        return None
    finally:
        if conn:
            conn.close()

def load_blast_species_set(filepath):
    with open(filepath, 'r') as f:
        return f.read().lower()

def compare_lca_and_blast(config_path, og_id, lca_files, blast_files, output_file):
    db_params = load_db_config(config_path)

    # Combine files
    concatenate_files(lca_files, f"lca_combined.{og_id}.tsv")
    concatenate_files(blast_files, f"blast_combined.{og_id}.tsv")

    # Get nominal_species_id from DB
    db_species = get_species_for_ogid(db_params, og_id)
    if db_species is None:
        print(f"[WARN] OG ID '{og_id}' nominal species not found in database.")
        return

    # Normalise once
    db_species_norm = normalise_name(db_species)

    # Load BLAST results
    blast_blob = load_blast_species_set(f"blast_combined.{og_id}.tsv")

    # Process LCA and compare
    with open(f"lca_combined.{og_id}.tsv", newline='') as tsvfile, open(output_file, "w", newline='') as out:
        reader = csv.reader(tsvfile, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["og_id", "LCA_result", "nom_species_id", "Match_YN", "Found_in_blast_YN"])

        for row in reader:
            if len(row) < 3:
                continue
            lca = row[2].strip()
            match = "Yes" if normalise_name(lca) == db_species_norm else "No"
            in_blast = "Yes" if db_species_norm in blast_blob else "No"
            writer.writerow([og_id, lca, db_species, match, in_blast])

    print(f"[INFO] Results written to: {output_file}")


# ---------------------------
# Entry point
# ---------------------------
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage:\n  python compare_lca_blast.py <db.cfg> <OG_ID> <lca1,lca2,...> <blast1,blast2,...>")
        sys.exit(1)

    config_file = sys.argv[1]
    og_id = sys.argv[2]
    lca_files = sys.argv[3].split(',')
    blast_files = sys.argv[4].split(',')

    output_file = f"lca_results.{og_id}.tsv"
    compare_lca_and_blast(config_file, og_id, lca_files, blast_files, output_file)
