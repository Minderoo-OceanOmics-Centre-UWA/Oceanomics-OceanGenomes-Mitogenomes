#!/usr/bin/env python3
#singularity run $SING/psycopg2:0.1.sif python push_lca_blast_results.py /home/tpeirce/postgresql_details/oceanomics.cfg OGtest lca_combined.tsv blast_combined.tsv
import psycopg2
import pandas as pd
import numpy as np
import argparse
import configparser
from pathlib import Path

def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    
    config = configparser.ConfigParser()
    config.read(config_file)

    required_keys = ['dbname', 'user', 'password', 'host', 'port']
    if not config.has_section('postgres'):
        raise ValueError("❌ Missing [postgres] section in config file.")
    for key in required_keys:
        if not config.has_option('postgres', key):
            raise ValueError(f"❌ Missing '{key}' in config file.")

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')
    }

blast_column_headers = [
    "query_id", "match_sequence_id", "taxon_id", "scientific_name",
    "common_name", "kingdoms", "percent_identity", "alignment_length",
    "query_length", "subject_length", "mismatch", "gap_open", "gaps",
    "query_start", "query_end", "subject_start", "subject_end",
    "subject_title", "evalue", "bit_score", "query_coverage", "subject_coverage",
    "blast_run_date", "region" 
]

def process_blast(blast_file, sample, db_params):
    print(f"📂 Reading BLAST file: {blast_file}")
    df = pd.read_csv(blast_file, sep='\t', header=None, names=blast_column_headers).replace({np.nan: None})

    if 'query_id' not in df.columns:
        print("❌ Missing 'query_id' in BLAST file")
        return

    # A sample whose regions produced no filtered hits gets an empty combined file. Splitting
    # query_id on an empty frame yields a 0-column result and raises "Columns must be same
    # length as key", which under errorStrategy 'ignore' silently dropped the sample's whole
    # BLAST upload. No hits is a legitimate outcome, not an error.
    if df.empty:
        print(f"ℹ️  No filtered BLAST hits for {sample} — nothing to upload")
        return

    df[['og_id', 'tech', 'seq_date', 'code', 'annotation']] = df['query_id'].str.split('.', expand=True)

    success, failure = 0, 0

    with psycopg2.connect(**db_params) as conn:
        cursor = conn.cursor()
        for _, row in df.iterrows():
            row_dict = row.to_dict()
            upsert_query = """
            INSERT INTO blast_filtered_lca (
                og_id, tech, seq_date, code, annotation, match_sequence_id, taxon_id,
                scientific_name, common_name, kingdoms, percent_identity, alignment_length,
                query_length, subject_length, mismatch, gap_open, gaps, query_start,
                query_end, subject_start, subject_end, subject_title, evalue, bit_score,
                query_coverage, subject_coverage, region, blast_run_date
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s, %(match_sequence_id)s, %(taxon_id)s,
                %(scientific_name)s, %(common_name)s, %(kingdoms)s, %(percent_identity)s, %(alignment_length)s,
                %(query_length)s, %(subject_length)s, %(mismatch)s, %(gap_open)s, %(gaps)s, %(query_start)s,
                %(query_end)s, %(subject_start)s, %(subject_end)s, %(subject_title)s, %(evalue)s, %(bit_score)s,
                %(query_coverage)s, %(subject_coverage)s, %(region)s, %(blast_run_date)s
            )
            ON CONFLICT (og_id, tech, seq_date, code, annotation, match_sequence_id, region)
            DO UPDATE SET
                blast_run_date = EXCLUDED.blast_run_date,
                taxon_id = EXCLUDED.taxon_id,
                scientific_name = EXCLUDED.scientific_name,
                common_name = EXCLUDED.common_name,
                kingdoms = EXCLUDED.kingdoms,
                percent_identity = EXCLUDED.percent_identity,
                alignment_length = EXCLUDED.alignment_length,
                query_length = EXCLUDED.query_length,
                subject_length = EXCLUDED.subject_length,
                mismatch = EXCLUDED.mismatch,
                gap_open = EXCLUDED.gap_open,
                gaps = EXCLUDED.gaps,
                query_start = EXCLUDED.query_start,
                query_end = EXCLUDED.query_end,
                subject_start = EXCLUDED.subject_start,
                subject_end = EXCLUDED.subject_end,
                subject_title = EXCLUDED.subject_title,
                evalue = EXCLUDED.evalue,
                bit_score = EXCLUDED.bit_score,
                query_coverage = EXCLUDED.query_coverage,
                subject_coverage = EXCLUDED.subject_coverage;
            """
            params = {
                "og_id": row_dict["og_id"],
                "tech": row_dict["tech"],
                "seq_date": row_dict["seq_date"],
                "code": row_dict["code"],
                "annotation": row_dict["annotation"],
                "match_sequence_id": row_dict.get("match_sequence_id"),
                "taxon_id": row_dict.get("taxon_id"),
                "scientific_name": row_dict.get("scientific_name"),
                "common_name": row_dict.get("common_name"),
                "kingdoms": row_dict.get("kingdoms"),
                "percent_identity": row_dict.get("percent_identity"),
                "alignment_length": row_dict.get("alignment_length"),
                "query_length": row_dict.get("query_length"),
                "subject_length": row_dict.get("subject_length"),
                "mismatch": row_dict.get("mismatch", 0),
                "gap_open": row_dict.get("gap_open", 0),
                "gaps": row_dict.get("gaps", 0),
                "query_start": row_dict.get("query_start"),
                "query_end": row_dict.get("query_end"),
                "subject_start": row_dict.get("subject_start"),
                "subject_end": row_dict.get("subject_end"),
                "subject_title": row_dict.get("subject_title"),
                "evalue": row_dict.get("evalue", 0),
                "bit_score": row_dict.get("bit_score"),
                "query_coverage": row_dict.get("query_coverage"),
                "subject_coverage": row_dict.get("subject_coverage"),
                "region": row_dict["region"],
                "blast_run_date": row_dict.get("blast_run_date")
            }

            try:
                cursor.execute(upsert_query, params)
                success += 1
            except Exception as e:
                failure += 1
                print(f"❌ BLAST row failed ({row_dict.get('query_id')}): {e}")

    print(f"✅ BLAST upload complete: {success} rows succeeded, {failure} failed")


lca_column_headers = [
    "seq_id", "species_in_LCA", "numberOfUnq_BlastHits", "domain",
    "phylum", "class", "order", "family", "genus", "specificEpithet",
    "scientificName", "scientificNameAuthorship", "taxonRank", "top_taxonID",
    "taxonID_db", "top_verbatimIdentification", "top_accession_id", "accession_id_ref_db",
    "taxonRank_db",
    "top_percent_match", "top_percent_query_cover", "top_percent_query_cover_hsp", "alignment_length",
    "subject_length", "sequence_length", "top_confidence_score", "sequence_region", "lca_run_date",
    "dna_sequence", "identificationRemarks"
]

# The group --force prunes: everything this sample/region has ever recorded.
LCA_HISTORY_KEY = ["og_id", "tech", "seq_date", "code", "annotation", "region"]


def prune_superseded_lca(cursor, written):
    """Delete lca rows this run did not produce, per sample/region it touched.

    `written` maps an LCA_HISTORY_KEY tuple to the set of lca_run_date values this
    run wrote for it. A row whose content matched has just had its date refreshed
    to the current run's, so anything left on another date is a superseded result.
    """
    where = " AND ".join(f'"{c}" = %s' for c in LCA_HISTORY_KEY)
    query = (
        f"DELETE FROM lca WHERE {where} "
        f'AND "lca_run_date"::text <> ALL(%s)'
    )
    deleted = 0
    for key, run_dates in written.items():
        cursor.execute(query, list(key) + [[str(d) for d in run_dates]])
        deleted += cursor.rowcount
    return deleted


def process_lca(lca_file, sample, db_params, force=False):
    print(f"📂 Reading LCA file: {lca_file}")
    # df = pd.read_csv(lca_file, sep='\t', header=True, names=lca_column_headers).replace({np.nan: None})
    # header=0 means a 0-byte file raises EmptyDataError rather than parsing to an empty
    # frame (process_blast passes explicit names, so it does not). Both mean the same thing
    # here -- nothing to upload -- so treat them the same instead of failing the task.
    try:
        df = pd.read_csv(lca_file, sep='\t', header=0).replace({np.nan: None})
    except pd.errors.EmptyDataError:
        print(f"ℹ️  Empty LCA file for {sample} — nothing to upload")
        return

    if 'seq_id' not in df.columns:
        print("❌ Missing 'seq_id' in LCA file")
        return

    # Same as process_blast: a header-only combined file is "no assignments", not a failure.
    if df.empty:
        print(f"ℹ️  No LCA assignments for {sample} — nothing to upload")
        return

    df[['og_id', 'tech', 'seq_date', 'code', 'annotation']] = df['seq_id'].str.split('.', expand=True)
    success, failure = 0, 0
    # HISTORY_KEY tuple -> set of lca_run_date values written for it, used by --force.
    written = {}

    with psycopg2.connect(**db_params) as conn:
        cursor = conn.cursor()
        for _, row in df.iterrows():
            row_dict = row.to_dict()
            upsert_query = """
            INSERT INTO lca (
                og_id,
                tech,
                seq_date,
                code,
                annotation,
                region,
                lca_run_date,
                species_in_lca,
                number_unq_blast_hits,
                domain,
                phylum,
                class,
                "order",
                family,
                genus,
                specific_epiphet,
                species,
                scientific_name_authorship,
                taxon_rank,
                top_taxon_id,
                taxon_id_db,
                top_accession_id,
                accession_id_ref_db,
                taxon_rank_db,
                top_percent_match,
                top_percent_query_cover,
                top_percent_query_cover_hsp,
                alignment_length,
                subject_length,
                sequence_length,
                top_confidence_score
            )
            VALUES (
                %(og_id)s,
                %(tech)s,
                %(seq_date)s,
                %(code)s,
                %(annotation)s,
                %(region)s,
                %(lca_run_date)s,
                %(species_in_lca)s,
                %(number_unq_blast_hits)s,
                %(domain)s,
                %(phylum)s,
                %(class)s,
                %(order)s,
                %(family)s,
                %(genus)s,
                %(specific_epiphet)s,
                %(species)s,
                %(scientific_name_authorship)s,
                %(taxon_rank)s,
                %(top_taxon_id)s,
                %(taxon_id_db)s,
                %(top_accession_id)s,
                %(accession_id_ref_db)s,
                %(taxon_rank_db)s,
                %(top_percent_match)s,
                %(top_percent_query_cover)s,
                %(top_percent_query_cover_hsp)s,
                %(alignment_length)s,
                %(subject_length)s,
                %(sequence_length)s,
                %(top_confidence_score)s
            )
            -- Rows are keyed by content (content_hash is set by the
            -- lca_content_hash trigger), so a conflict means this run reproduced
            -- an existing result exactly. The only thing left to record is that
            -- it was seen again; a differing result hashes differently and is
            -- inserted as a new row rather than overwriting the old one.
            ON CONFLICT ON CONSTRAINT lca_content_unique
            DO UPDATE SET
                lca_run_date = EXCLUDED.lca_run_date;
            """

            # params: 1:1 mapping between placeholders and df column names
            params = {
                "og_id": row_dict.get("og_id"),
                "tech": row_dict.get("tech"),
                "seq_date": row_dict.get("seq_date"),
                "code": row_dict.get("code"),
                "annotation": row_dict.get("annotation"),
                "region": row_dict.get("sequence_region"),
                "lca_run_date": row_dict.get("lca_run_date"),

                "species_in_lca": row_dict.get("species_in_LCA"),
                "number_unq_blast_hits": row_dict.get("numberOfUnq_BlastHits"),
                "domain": row_dict.get("domain"),
                "phylum": row_dict.get("phylum"),
                "class": row_dict.get("class"),
                "order": row_dict.get("order"),
                "family": row_dict.get("family"),
                "genus": row_dict.get("genus"),

                "specific_epiphet": row_dict.get("specificEpithet"),
                "species": row_dict.get("scientificName"),
                "scientific_name_authorship": row_dict.get("scientificNameAuthorship"),
                "taxon_rank": row_dict.get("taxonRank"),

                "top_taxon_id": row_dict.get("top_taxonID"),
                "taxon_id_db": row_dict.get("taxonID_db"),
                "top_accession_id": row_dict.get("top_accession_id"),
                "accession_id_ref_db": row_dict.get("accession_id_ref_db"),
                "taxon_rank_db": row_dict.get("taxonRank_db"),
                "top_percent_match": row_dict.get("top_percent_match"),
                "top_percent_query_cover": row_dict.get("top_percent_query_cover"),
                "top_percent_query_cover_hsp": row_dict.get("top_percent_query_cover_hsp"),

                "alignment_length": row_dict.get("alignment_length"),
                "subject_length": row_dict.get("subject_length"),
                "sequence_length": row_dict.get("sequence_length"),
                "top_confidence_score": row_dict.get("top_confidence_score"),
            }

            try:
                cursor.execute(upsert_query, params)
                success += 1
                key = tuple(params[c] for c in LCA_HISTORY_KEY)
                written.setdefault(key, set()).add(params["lca_run_date"])
            except Exception as e:
                failure += 1
                print(f"❌ LCA row failed ({row_dict.get('seq_id')}): {e}")

        # Prune only after every row landed: deleting around a failure could drop
        # the older result without its replacement having been written.
        if force and failure == 0:
            pruned = prune_superseded_lca(cursor, written)
            print(f"🧹 --force: deleted {pruned} superseded lca rows")
        elif force:
            print("⚠️ Skipping --force prune: upload had failures.")

    print(f"✅ LCA upload complete: {success} rows succeeded, {failure} failed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Upload filtered BLAST hits and LCA assignments for one sample."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "After uploading, delete lca rows for the same sample/region that this "
            "run did not produce, so only the current result is kept. Off by "
            "default: superseded results are preserved alongside the new ones. "
            "Does not affect blast_filtered_lca, which overwrites in place and "
            "keeps no history."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("sample")
    parser.add_argument("lca_results")
    parser.add_argument("blast_results")
    args = parser.parse_args()

    db_config = load_db_config(args.config_file)

    print("\n🔄 Starting BLAST processing...")
    process_blast(args.blast_results, args.sample, db_config)

    print("\n🔄 Starting LCA processing...")
    process_lca(args.lca_results, args.sample, db_config, force=args.force)

    print("\n🎉 All processing complete.")
