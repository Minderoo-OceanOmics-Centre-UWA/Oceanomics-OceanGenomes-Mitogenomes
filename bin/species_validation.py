#!/usr/bin/env python3
import argparse
import psycopg2
import csv
import configparser
import sys

SPECIES_IN_LCA_COLUMN = "species_in_LCA"   # <--- NEW

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

def concatenate_lca_files(file_list, output_file):
    """
    Concatenate LCA files that each have a header line.
    Writes the header only once (from the first file).
    """
    first_file = True
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                for line_num, line in enumerate(infile):
                    # Always write the header from the first file
                    if first_file:
                        outfile.write(line)
                    else:
                        # Skip header line for subsequent files
                        if line_num == 0:
                            continue
                        outfile.write(line)
            first_file = False

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
        # We keep this as a big lowercased blob and search for the
        # normalised DB species name in it later.
        return f.read().lower()

def parse_assembly_key_from_blast(blast_file):
    """
    Extract (og_id, tech, seq_date, code, annotation) from the first
    query_id field in a BLAST table. The query_id is expected to look like
    'og_id.tech.seq_date.code.annotation' (5 dot-separated parts).
    Returns (og_id, tech, seq_date, code, annotation) or None if the file
    is empty / malformed.
    """
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                first_field = line.split('\t', 1)[0]
                parts = first_field.split('.')
                if len(parts) >= 5:
                    return tuple(parts[:5])
                return None
    except Exception as e:
        print(f"[WARN] Could not parse assembly key from BLAST file '{blast_file}': {e}")
    return None


def upsert_lca_validation(
    db_params, key, validated_species_name, validator="nf-core", force=False
):
    """
    Insert/update the lca_validation row for this sample with the validated
    species name and validator tag. The composite key matches the assembly
    naming convention: (og_id, tech, seq_date, code, annotation).

    Guard: if an existing row carries a validator tag other than "nf-core",
    the row was probably set by a human reviewer or another pipeline and
    must not be silently overwritten by an automated re-run. Skip the
    upsert in that case unless ``force=True``.
    """
    og_id, tech, seq_date, code, annotation = key
    params = {
        "og_id": og_id,
        "tech": tech,
        "seq_date": seq_date,
        "code": code,
        "annotation": annotation,
        "validated_species_name": validated_species_name,
        "validator": validator,
    }
    conn = None
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            if not force:
                cur.execute(
                    """
                    SELECT validator, validated_species_name
                    FROM lca_validation
                    WHERE og_id = %(og_id)s
                      AND tech = %(tech)s
                      AND seq_date = %(seq_date)s
                      AND code = %(code)s
                      AND annotation = %(annotation)s
                    """,
                    params,
                )
                existing = cur.fetchone()
                if existing is not None:
                    existing_validator, existing_species = existing
                    if (
                        existing_validator is not None
                        and str(existing_validator).strip().lower() != "nf-core"
                    ):
                        print(
                            "⚠️ Existing values preserved for lca_validation "
                            f"{og_id}.{tech}.{seq_date}.{code}.{annotation}: "
                            f"validator='{existing_validator}' (not 'nf-core'); "
                            "pass --force to overwrite."
                        )
                        print(
                            "📌 Preserved stored values: "
                            f"{{'validator': '{existing_validator}', "
                            f"'validated_species_name': '{existing_species}'}}"
                        )
                        return True

            upsert_query = """
            INSERT INTO lca_validation (
                og_id, tech, seq_date, code, annotation,
                validated_species_name, validator
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s,
                %(validated_species_name)s, %(validator)s
            )
            ON CONFLICT (og_id, tech, seq_date, code, annotation)
            DO UPDATE SET
                validated_species_name = EXCLUDED.validated_species_name,
                validator = EXCLUDED.validator
            """
            cur.execute(upsert_query, params)
        conn.commit()
        if force:
            print(
                "✅ Success: lca_validation overwritten (--force) for "
                f"{og_id}.{tech}.{seq_date}.{code}.{annotation} "
                f"-> validated_species_name='{validated_species_name}', validator='{validator}'"
            )
        else:
            print(
                "✅ Success: lca_validation upserted for "
                f"{og_id}.{tech}.{seq_date}.{code}.{annotation} "
                f"-> validated_species_name='{validated_species_name}', validator='{validator}'"
            )
        return True
    except Exception as e:
        if conn is not None:
            conn.rollback()
        print(f"❌ Database error while writing lca_validation: {e}")
        return False
    finally:
        if conn is not None:
            conn.close()


def compare_lca_and_blast(config_path, og_id, lca_files, blast_files, output_file, assembly_prefix=None, force=False):
    db_params = load_db_config(config_path)
    # File-naming prefix only: an OG can have multiple assembly attempts, so
    # combined/summary filenames must be qualified with the unique per-assembly
    # prefix (meta.mt_assembly_prefix) to avoid colliding with other attempts
    # for the same OG. The DB lookups and output row content still use the
    # real og_id.
    prefix = assembly_prefix or og_id

    # Combine files
    concatenate_lca_files(lca_files, f"lca_combined.{prefix}.tsv")
    concatenate_files(blast_files, f"blast_combined.{prefix}.tsv")

    # Get nominal_species_id from DB
    db_species = get_species_for_ogid(db_params, og_id)
    if db_species is None:
        print(f"[WARN] OG ID '{og_id}' nominal species not found in database.")
        return

    # Normalise once
    db_species_norm = normalise_name(db_species)

    # Load BLAST results
    blast_blob = load_blast_species_set(f"blast_combined.{prefix}.tsv")

    sample_validated = False

    # Process LCA and compare
    with open(f"lca_combined.{prefix}.tsv", newline='') as tsvfile, open(output_file, "w", newline='') as out:
        # Use DictReader so we can refer to columns by name
        reader = csv.DictReader(tsvfile, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["og_id", "LCA_result", "nom_species_id", "Match_YN", "Found_in_blast_YN"])

        for row in reader:
            if not row:
                continue

            # Get the comma-separated species list from the LCA file
            species_in_lca_raw = row.get(SPECIES_IN_LCA_COLUMN)
            if not species_in_lca_raw:
                # If the column isn't present or is empty, skip this row
                continue

            # Split on commas and normalise each species name
            species_list_norm = [
                normalise_name(s)
                for s in species_in_lca_raw.split(',')
                if s.strip()
            ]

            # Check if DB nominal species is among the LCA species list
            match = "Yes" if db_species_norm in species_list_norm else "No"

            # Check if nominal species appears in BLAST output (as before)
            in_blast = "Yes" if db_species_norm in blast_blob else "No"
            if in_blast == "Yes":
                sample_validated = True

            # LCA_result column: store the raw species_in_LCA string
            writer.writerow([og_id, species_in_lca_raw, db_species, match, in_blast])

    print(f"[INFO] Results written to: {output_file}")

    # Push the validation result to the lca_validation table. We only write a
    # row when the sample is validated (Found_in_blast_YN = Yes for at least
    # one region) so unvalidated samples leave the existing DB row untouched.
    if sample_validated:
        key = parse_assembly_key_from_blast(f"blast_combined.{prefix}.tsv")
        if key is None:
            print(
                "❌ Could not derive (og_id, tech, seq_date, code, annotation) from "
                f"blast_combined.{prefix}.tsv — skipping lca_validation upsert."
            )
        else:
            upsert_lca_validation(db_params, key, db_species, validator="nf-core", force=force)
    else:
        print(
            f"ℹ️ Sample {og_id} not validated (no Found_in_blast_YN=Yes) — "
            "skipping lca_validation upsert."
        )


# ---------------------------
# Entry point
# ---------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compare per-region LCA / BLAST results against the nominal species "
            "stored in the OceanOmics DB, write a summary TSV, and (when "
            "validated) upsert the result into the lca_validation table."
        )
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite an existing lca_validation row even when its validator "
            "is something other than 'nf-core' (e.g. a human reviewer). Off "
            "by default."
        ),
    )
    parser.add_argument(
        "--assembly-prefix",
        default=None,
        help=(
            "Unique per-assembly filename prefix (meta.mt_assembly_prefix). "
            "Used only to name the combined/summary output files so multiple "
            "assembly attempts for the same OG don't collide; defaults to "
            "og_id when omitted."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("og_id")
    parser.add_argument(
        "lca_files",
        help="Comma-separated list of per-region LCA TSVs.",
    )
    parser.add_argument(
        "blast_files",
        help="Comma-separated list of filtered BLAST TSVs.",
    )
    args = parser.parse_args()

    lca_files = args.lca_files.split(',')
    blast_files = args.blast_files.split(',')

    prefix = args.assembly_prefix or args.og_id
    output_file = f"lca_results.{prefix}.tsv"
    compare_lca_and_blast(
        args.config_file,
        args.og_id,
        lca_files,
        blast_files,
        output_file,
        assembly_prefix=args.assembly_prefix,
        force=args.force,
    )
