#!/usr/bin/env python3
#singularity run $SING/psycopg2:0.1.sif python push_emma_annotation_results.py /home/tpeirce/postgresql_details/oceanomics.cfg OGtest /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mitogenomes/OG900/OG900.ilmn.250131.getorg1770/annotation/structure_check/annotation_stats.csv
import argparse
import psycopg2
import pandas as pd
import numpy as np
import configparser
import sys
import re
from pathlib import Path

# -------------------------------
# Load DB credentials from .cfg
# -------------------------------
def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    
    config = configparser.ConfigParser()
    config.read(config_file)

    if not config.has_section('postgres'):
        raise ValueError("❌ Missing [postgres] section in config file.")

    required_keys = ['dbname', 'user', 'password', 'host', 'port']
    for key in required_keys:
        if not config.has_option('postgres', key):
            raise ValueError(f"❌ Missing '{key}' in [postgres] section of config file.")

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')    
    }

# -------------------------------
# Main logic
# -------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Push EMMA annotation stats to the OceanOmics DB."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite the annotation columns even when the existing row "
            "already has some of them filled. Off by default: a row is only "
            "(re)written when every annotation cell is empty; partially-filled "
            "rows are preserved."
        ),
    )
    parser.add_argument("config_file")
    parser.add_argument("sample")
    parser.add_argument("annotation_stats")
    args = parser.parse_args()

    config_file = args.config_file
    sample = args.sample
    annotation_stats = args.annotation_stats
    force_overwrite = args.force

    print(f"Importing data from {annotation_stats}")
    # Load and preprocess annotation data
    stats = pd.read_csv(annotation_stats)

    # Replace NaN with None
    stats = stats.replace({np.nan: None})

    # Normalize column names (remove spaces & make lowercase)
    stats.columns = stats.columns.str.strip().str.lower()

    try:
        db_params = load_db_config(config_file)
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()

        # Annotation columns owned by this script. The PK row is pre-created by
        # the assembly upload, so we don't gate on row existence — we gate on
        # whether these cells are already (partially) filled.
        annotation_columns = [
            "annotation", "extra_genes", "missing_genes", "order_correct", "passed", "length_emma", "seqlength_12s",
            "seqlength_16s", "seqlength_co1", "cds_no", "trna_no", "rrna_no", "rrna12s",
            "rrna16s", "atp6", "atp8", "cox1", "cox2", "cox3", "cytb", "nad1", "nad2", "nad3", "nad4", "nad4l",
            "nad5", "nad6", "trna_phe", "trna_val", "trna_leuuag", "trna_leuuaa", "trna_ile", "trna_met",
            "trna_thr", "trna_pro", "trna_lys", "trna_asp", "trna_glu", "trna_sergcu", "trna_seruga",
            "trna_tyr", "trna_cys", "trna_trp", "trna_ala", "trna_asn", "trna_gly", "trna_arg", "trna_his",
            "trna_gln", "atp6_trans", "atp8_trans", "cox1_trans", "cox2_trans", "cox3_trans", "cytb_trans", "nad1_trans",
            "nad2_trans", "nad3_trans", "nad4_trans", "nad4l_trans", "nad5_trans", "nad6_trans"
        ]

        for index, row in stats.iterrows():
            row_dict = row.to_dict()
            # Extract primary key values
            og_id, tech, seq_date, code = row['og_id'], row['tech'], row['seq_date'], row['code']
            # seq_date is a TEXT column in the DB, but pandas reads 260514 as an int.
            # Bind it as a string so the existence-check comparison (text = %s) doesn't
            # fail with "operator does not exist: text = integer".
            if seq_date is not None:
                seq_date = str(seq_date)

            upsert_query = """
            INSERT INTO mitogenome_data (
                og_id, tech, seq_date, code, annotation, extra_genes, missing_genes, order_correct, passed, length_emma, seqlength_12s,
                seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, rrna12s,
                rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, 
                nad5, nad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, 
                trna_thr, trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, 
                trna_tyr, trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, 
                trna_gln, atp6_trans, atp8_trans, cox1_trans, cox2_trans, cox3_trans, cytb_trans, nad1_trans, 
                nad2_trans, nad3_trans, nad4_trans, nad4l_trans, nad5_trans, nad6_trans
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s, %(extra_genes)s, %(missing_genes)s, %(order_correct)s, %(passed)s, %(length_emma)s, %(seqlength_12s)s,
                %(seqlength_16s)s, %(seqlength_co1)s, %(cds_no)s, %(trna_no)s, %(rrna_no)s, %(rrna12s)s,
                %(rrna16s)s, %(atp6)s, %(atp8)s, %(cox1)s, %(cox2)s, %(cox3)s, %(cytb)s, %(nad1)s, %(nad2)s, %(nad3)s, %(nad4)s, %(nad4l)s,
                %(nad5)s, %(nad6)s, %(trna_phe)s, %(trna_val)s, %(trna_leuuag)s, %(trna_leuuaa)s, %(trna_ile)s, %(trna_met)s,
                %(trna_thr)s, %(trna_pro)s, %(trna_lys)s, %(trna_asp)s, %(trna_glu)s, %(trna_sergcu)s, %(trna_seruga)s,
                %(trna_tyr)s, %(trna_cys)s, %(trna_trp)s, %(trna_ala)s, %(trna_asn)s, %(trna_gly)s, %(trna_arg)s, %(trna_his)s,
                %(trna_gln)s, %(atp6_trans)s, %(atp8_trans)s, %(cox1_trans)s, %(cox2_trans)s, %(cox3_trans)s, %(cytb_trans)s, %(nad1_trans)s, %(nad2_trans)s, 
                %(nad3_trans)s, %(nad4_trans)s, %(nad4l_trans)s, %(nad5_trans)s, %(nad6_trans)s
            )
            ON CONFLICT (og_id, tech, seq_date, code) 
            DO UPDATE SET
                annotation = EXCLUDED.annotation,
                extra_genes = EXCLUDED.extra_genes,
                missing_genes = EXCLUDED.missing_genes,
                order_correct = EXCLUDED.order_correct,
                passed = EXCLUDED.passed,
                length_emma = EXCLUDED.length_emma,
                seqlength_12s = EXCLUDED.seqlength_12s,
                seqlength_16s = EXCLUDED.seqlength_16s,
                seqlength_co1 = EXCLUDED.seqlength_co1,
                cds_no = EXCLUDED.cds_no,
                trna_no = EXCLUDED.trna_no,
                rrna_no = EXCLUDED.rrna_no,
                rrna12s = EXCLUDED.rrna12s,
                rrna16s = EXCLUDED.rrna16s,
                atp6 = EXCLUDED.atp6,
                atp8 = EXCLUDED.atp8,
                cox1 = EXCLUDED.cox1,
                cox2 = EXCLUDED.cox2,
                cox3 = EXCLUDED.cox3,
                cytb = EXCLUDED.cytb,
                nad1 = EXCLUDED.nad1,
                nad2 = EXCLUDED.nad2,
                nad3 = EXCLUDED.nad3,
                nad4 = EXCLUDED.nad4,
                nad4l = EXCLUDED.nad4l,
                nad5 = EXCLUDED.nad5,
                nad6 = EXCLUDED.nad6,
                trna_phe = EXCLUDED.trna_phe,
                trna_val = EXCLUDED.trna_val,
                trna_leuuag = EXCLUDED.trna_leuuag,
                trna_leuuaa = EXCLUDED.trna_leuuaa,
                trna_ile = EXCLUDED.trna_ile,
                trna_met = EXCLUDED.trna_met,
                trna_thr = EXCLUDED.trna_thr,
                trna_pro = EXCLUDED.trna_pro,
                trna_lys = EXCLUDED.trna_lys,
                trna_asp = EXCLUDED.trna_asp,
                trna_glu = EXCLUDED.trna_glu,
                trna_sergcu = EXCLUDED.trna_sergcu,
                trna_seruga = EXCLUDED.trna_seruga,
                trna_tyr = EXCLUDED.trna_tyr,
                trna_cys = EXCLUDED.trna_cys,
                trna_trp = EXCLUDED.trna_trp,
                trna_ala = EXCLUDED.trna_ala,
                trna_asn = EXCLUDED.trna_asn,
                trna_gly = EXCLUDED.trna_gly,
                trna_arg = EXCLUDED.trna_arg,
                trna_his = EXCLUDED.trna_his,
                trna_gln = EXCLUDED.trna_gln,
                atp6_trans = EXCLUDED.atp6_trans,
                atp8_trans = EXCLUDED.atp8_trans,
                cox1_trans = EXCLUDED.cox1_trans,
                cox2_trans = EXCLUDED.cox2_trans,
                cox3_trans = EXCLUDED.cox3_trans,
                cytb_trans = EXCLUDED.cytb_trans,
                nad1_trans = EXCLUDED.nad1_trans,
                nad2_trans = EXCLUDED.nad2_trans,
                nad3_trans = EXCLUDED.nad3_trans,
                nad4_trans = EXCLUDED.nad4_trans,
                nad4l_trans = EXCLUDED.nad4l_trans,
                nad5_trans = EXCLUDED.nad5_trans,
                nad6_trans = EXCLUDED.nad6_trans
            RETURNING annotation, extra_genes, missing_genes, order_correct, passed, length_emma, seqlength_12s,
                seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, rrna12s,
                rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, 
                nad5, nad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, 
                trna_thr, trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, 
                trna_tyr, trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, 
                trna_gln, atp6_trans, atp8_trans, cox1_trans, cox2_trans, cox3_trans, cytb_trans, nad1_trans, 
                nad2_trans, nad3_trans, nad4_trans, nad4l_trans, nad5_trans, nad6_trans
            """

            params = {
                "og_id": row_dict["og_id"],
                "tech": row_dict["tech"],
                "seq_date": seq_date,
                "code": row_dict["code"],
                "annotation": row_dict["annotation"],
                "extra_genes": row_dict["extra_genes"],
                "missing_genes": row_dict["missing_genes"],
                "order_correct": row_dict["order_correct"],
                "passed": row_dict["passed"],
                "length_emma": row_dict.get("total_length"),
                "seqlength_12s": row_dict.get("seqlength_12s"),
                "seqlength_16s": row_dict.get("seqlength_16s"),
                "seqlength_co1": row_dict.get("seqlength_co1"),
                "cds_no": row_dict.get("num_cds"),
                "trna_no": row_dict.get("num_trna"),
                "rrna_no": row_dict.get("num_rrna"),
                "rrna12s": row_dict.get("rnr1"),
                "rrna16s": row_dict.get("rnr2"),
                "atp6": row_dict.get("atp6"),
                "atp8": row_dict.get("atp8"),
                "cox1": row_dict.get("co1"),
                "cox2": row_dict.get("co2"),
                "cox3": row_dict.get("co3"),
                "cytb": row_dict.get("cytb"),
                "nad1": row_dict.get("nd1"),
                "nad2": row_dict.get("nd2"),
                "nad3": row_dict.get("nd3"),
                "nad4": row_dict.get("nd4"),
                "nad4l": row_dict.get("nd4l"),
                "nad5": row_dict.get("nd5"),
                "nad6": row_dict.get("nd6"),
                "trna_phe": row_dict.get("tf"),
                "trna_val": row_dict.get("tv"),
                "trna_leuuag": row_dict.get("tl1"),
                "trna_leuuaa": row_dict.get("tl2"),
                "trna_ile": row_dict.get("ti"),
                "trna_met": row_dict.get("tm"),
                "trna_thr": row_dict.get("tt"),
                "trna_pro": row_dict.get("tp"),
                "trna_lys": row_dict.get("tk"),
                "trna_asp": row_dict.get("td"),
                "trna_glu": row_dict.get("te"),
                "trna_sergcu": row_dict.get("ts1"),
                "trna_seruga": row_dict.get("ts2"),
                "trna_tyr": row_dict.get("ty"),
                "trna_cys": row_dict.get("tc"),
                "trna_trp": row_dict.get("tw"),
                "trna_ala": row_dict.get("ta"),
                "trna_asn": row_dict.get("tn"),
                "trna_gly": row_dict.get("tg"),
                "trna_arg": row_dict.get("tr"),
                "trna_his": row_dict.get("th"),
                "trna_gln": row_dict.get("tq"),
                "atp6_trans": row_dict.get("atp6_trans"),
                "atp8_trans": row_dict.get("atp8_trans"),
                "cox1_trans": row_dict.get("co1_trans"),
                "cox2_trans": row_dict.get("co2_trans"),
                "cox3_trans": row_dict.get("co3_trans"),
                "cytb_trans": row_dict.get("cytb_trans"),
                "nad1_trans": row_dict.get("nd1_trans"),
                "nad2_trans": row_dict.get("nd2_trans"),
                "nad3_trans": row_dict.get("nd3_trans"),
                "nad4_trans": row_dict.get("nd4_trans"),
                "nad4l_trans": row_dict.get("nd4l_trans"),
                "nad5_trans": row_dict.get("nd5_trans"),
                "nad6_trans": row_dict.get("nd6_trans"),
            }

            # Only (re)write when every annotation cell is empty. If any cell is
            # already filled (e.g. a partial annotation from a previous run), we
            # preserve the whole row so we never end up with a mix of old and new
            # values. --force bypasses this and overwrites unconditionally.
            if not force_overwrite:
                cursor.execute(
                    f"""
                    SELECT {", ".join(annotation_columns)}
                    FROM mitogenome_data
                    WHERE og_id = %(og_id)s
                      AND tech = %(tech)s
                      AND seq_date = %(seq_date)s
                      AND code = %(code)s
                    """,
                    {"og_id": og_id, "tech": tech, "seq_date": seq_date, "code": code},
                )
                existing = cursor.fetchone()
                if existing is not None and any(v is not None for v in existing):
                    preserved = dict(zip(annotation_columns, existing))
                    attempted = {col: params[col] for col in annotation_columns}
                    print(
                        f"⚠️ Existing values preserved for {og_id}.{tech}.{seq_date}.{code}: "
                        "row already has annotation data; pass --force to overwrite."
                    )
                    print(f"   🔒 Preserved (kept) values:    {preserved}")
                    print(f"   ⬆️ New values NOT uploaded:    {attempted}")
                    continue

            cursor.execute(upsert_query, params)
            returned = cursor.fetchone()
            conn.commit()

            print(f"📌 Final stored values: {dict(zip(annotation_columns, returned))}")
            if force_overwrite:
                print(f"✅ Success: Overwrote mitogenome_data for {sample} (--force).")
            else:
                print(f"✅ Success: Inserted/Updated mitogenome_data for {sample}")


    except Exception as e:
        if 'conn' in locals():
            conn.rollback()
        print(f"❌ Database error: {e}")

    finally:
        if 'cursor' in locals():
            cursor.close()
        if 'conn' in locals():
            conn.close()
