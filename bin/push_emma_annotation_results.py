#!/usr/bin/env python3
#singularity run $SING/psycopg2:0.1.sif python push_emma_annotation_results.py /home/tpeirce/postgresql_details/oceanomics.cfg OGtest /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR/mitogenomes/OG900/OG900.ilmn.250131.getorg1770/emma/structure_check/annotation_stats.csv
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
    if len(sys.argv) != 4:
        print("Usage:\n  push_emma_annotation_results.py <config_file> <meta.id> <annotation_stats>")
        sys.exit(1)

    config_file = sys.argv[1]
    sample = sys.argv[2]
    annotation_stats = sys.argv[3]

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
        
        for index, row in stats.iterrows():
            row_dict = row.to_dict()
            # Extract primary key values
            og_id, tech, seq_date, code = row['og_id'], row['tech'], row['seq_date'], row['code']

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
                annotation = CASE WHEN mitogenome_data.annotation IS NULL THEN EXCLUDED.annotation ELSE mitogenome_data.annotation END,
                extra_genes = CASE WHEN mitogenome_data.extra_genes IS NULL THEN EXCLUDED.extra_genes ELSE mitogenome_data.extra_genes END,
                missing_genes = CASE WHEN mitogenome_data.missing_genes IS NULL THEN EXCLUDED.missing_genes ELSE mitogenome_data.missing_genes END,
                order_correct = CASE WHEN mitogenome_data.order_correct IS NULL THEN EXCLUDED.order_correct ELSE mitogenome_data.order_correct END,
                passed = CASE WHEN mitogenome_data.passed IS NULL THEN EXCLUDED.passed ELSE mitogenome_data.passed END,
                length_emma = CASE WHEN mitogenome_data.length_emma IS NULL THEN EXCLUDED.length_emma ELSE mitogenome_data.length_emma END,
                seqlength_12s = CASE WHEN mitogenome_data.seqlength_12s IS NULL THEN EXCLUDED.seqlength_12s ELSE mitogenome_data.seqlength_12s END,
                seqlength_16s = CASE WHEN mitogenome_data.seqlength_16s IS NULL THEN EXCLUDED.seqlength_16s ELSE mitogenome_data.seqlength_16s END,
                seqlength_co1 = CASE WHEN mitogenome_data.seqlength_co1 IS NULL THEN EXCLUDED.seqlength_co1 ELSE mitogenome_data.seqlength_co1 END,
                cds_no = CASE WHEN mitogenome_data.cds_no IS NULL THEN EXCLUDED.cds_no ELSE mitogenome_data.cds_no END,
                trna_no = CASE WHEN mitogenome_data.trna_no IS NULL THEN EXCLUDED.trna_no ELSE mitogenome_data.trna_no END,
                rrna_no = CASE WHEN mitogenome_data.rrna_no IS NULL THEN EXCLUDED.rrna_no ELSE mitogenome_data.rrna_no END,
                rrna12s = CASE WHEN mitogenome_data.rrna12s IS NULL THEN EXCLUDED.rrna12s ELSE mitogenome_data.rrna12s END,
                rrna16s = CASE WHEN mitogenome_data.rrna16s IS NULL THEN EXCLUDED.rrna16s ELSE mitogenome_data.rrna16s END,
                atp6  = CASE WHEN mitogenome_data.atp6 IS NULL THEN EXCLUDED.atp6 ELSE mitogenome_data.atp6 END,
                atp8  = CASE WHEN mitogenome_data.atp8 IS NULL THEN EXCLUDED.atp8 ELSE mitogenome_data.atp8 END,
                cox1  = CASE WHEN mitogenome_data.cox1 IS NULL THEN EXCLUDED.cox1 ELSE mitogenome_data.cox1 END,
                cox2  = CASE WHEN mitogenome_data.cox2 IS NULL THEN EXCLUDED.cox2 ELSE mitogenome_data.cox2 END,
                cox3  = CASE WHEN mitogenome_data.cox3 IS NULL THEN EXCLUDED.cox3 ELSE mitogenome_data.cox3 END,
                cytb  = CASE WHEN mitogenome_data.cytb IS NULL THEN EXCLUDED.cytb ELSE mitogenome_data.cytb END,
                nad1  = CASE WHEN mitogenome_data.nad1 IS NULL THEN EXCLUDED.nad1 ELSE mitogenome_data.nad1 END,
                nad2  = CASE WHEN mitogenome_data.nad2 IS NULL THEN EXCLUDED.nad2 ELSE mitogenome_data.nad2 END,
                nad3  = CASE WHEN mitogenome_data.nad3 IS NULL THEN EXCLUDED.nad3 ELSE mitogenome_data.nad3 END,
                nad4  = CASE WHEN mitogenome_data.nad4 IS NULL THEN EXCLUDED.nad4 ELSE mitogenome_data.nad4 END,
                nad4l = CASE WHEN mitogenome_data.nad4l IS NULL THEN EXCLUDED.nad4l ELSE mitogenome_data.nad4l END,
                nad5  = CASE WHEN mitogenome_data.nad5 IS NULL THEN EXCLUDED.nad5 ELSE mitogenome_data.nad5 END,
                nad6  = CASE WHEN mitogenome_data.nad6 IS NULL THEN EXCLUDED.nad6 ELSE mitogenome_data.nad6 END,
                trna_phe = CASE WHEN mitogenome_data.trna_phe IS NULL THEN EXCLUDED.trna_phe ELSE mitogenome_data.trna_phe END,
                trna_val = CASE WHEN mitogenome_data.trna_val IS NULL THEN EXCLUDED.trna_val ELSE mitogenome_data.trna_val END,
                trna_leuuag = CASE WHEN mitogenome_data.trna_leuuag IS NULL THEN EXCLUDED.trna_leuuag ELSE mitogenome_data.trna_leuuag END,
                trna_leuuaa = CASE WHEN mitogenome_data.trna_leuuaa IS NULL THEN EXCLUDED.trna_leuuaa ELSE mitogenome_data.trna_leuuaa END,
                trna_ile = CASE WHEN mitogenome_data.trna_ile IS NULL THEN EXCLUDED.trna_ile ELSE mitogenome_data.trna_ile END,
                trna_met = CASE WHEN mitogenome_data.trna_met IS NULL THEN EXCLUDED.trna_met ELSE mitogenome_data.trna_met END,
                trna_thr = CASE WHEN mitogenome_data.trna_thr IS NULL THEN EXCLUDED.trna_thr ELSE mitogenome_data.trna_thr END,
                trna_pro = CASE WHEN mitogenome_data.trna_pro IS NULL THEN EXCLUDED.trna_pro ELSE mitogenome_data.trna_pro END,
                trna_lys = CASE WHEN mitogenome_data.trna_lys IS NULL THEN EXCLUDED.trna_lys ELSE mitogenome_data.trna_lys END,
                trna_asp = CASE WHEN mitogenome_data.trna_asp IS NULL THEN EXCLUDED.trna_asp ELSE mitogenome_data.trna_asp END,
                trna_glu = CASE WHEN mitogenome_data.trna_glu IS NULL THEN EXCLUDED.trna_glu ELSE mitogenome_data.trna_glu END,
                trna_sergcu = CASE WHEN mitogenome_data.trna_sergcu IS NULL THEN EXCLUDED.trna_sergcu ELSE mitogenome_data.trna_sergcu END,
                trna_seruga = CASE WHEN mitogenome_data.trna_seruga IS NULL THEN EXCLUDED.trna_seruga ELSE mitogenome_data.trna_seruga END,
                trna_tyr = CASE WHEN mitogenome_data.trna_tyr IS NULL THEN EXCLUDED.trna_tyr ELSE mitogenome_data.trna_tyr END,
                trna_cys = CASE WHEN mitogenome_data.trna_cys IS NULL THEN EXCLUDED.trna_cys ELSE mitogenome_data.trna_cys END,
                trna_trp = CASE WHEN mitogenome_data.trna_trp IS NULL THEN EXCLUDED.trna_trp ELSE mitogenome_data.trna_trp END,
                trna_ala = CASE WHEN mitogenome_data.trna_ala IS NULL THEN EXCLUDED.trna_ala ELSE mitogenome_data.trna_ala END,
                trna_asn = CASE WHEN mitogenome_data.trna_asn IS NULL THEN EXCLUDED.trna_asn ELSE mitogenome_data.trna_asn END,
                trna_gly = CASE WHEN mitogenome_data.trna_gly IS NULL THEN EXCLUDED.trna_gly ELSE mitogenome_data.trna_gly END,
                trna_arg = CASE WHEN mitogenome_data.trna_arg IS NULL THEN EXCLUDED.trna_arg ELSE mitogenome_data.trna_arg END,
                trna_his = CASE WHEN mitogenome_data.trna_his IS NULL THEN EXCLUDED.trna_his ELSE mitogenome_data.trna_his END,
                trna_gln = CASE WHEN mitogenome_data.trna_gln IS NULL THEN EXCLUDED.trna_gln ELSE mitogenome_data.trna_gln END,
                atp6_trans = CASE WHEN mitogenome_data.atp6_trans IS NULL THEN EXCLUDED.atp6_trans ELSE mitogenome_data.atp6_trans END,
                atp8_trans = CASE WHEN mitogenome_data.atp8_trans IS NULL THEN EXCLUDED.atp8_trans ELSE mitogenome_data.atp8_trans END,
                cox1_trans = CASE WHEN mitogenome_data.cox1_trans IS NULL THEN EXCLUDED.cox1_trans ELSE mitogenome_data.cox1_trans END,
                cox2_trans = CASE WHEN mitogenome_data.cox2_trans IS NULL THEN EXCLUDED.cox2_trans ELSE mitogenome_data.cox2_trans END,
                cox3_trans = CASE WHEN mitogenome_data.cox3_trans IS NULL THEN EXCLUDED.cox3_trans ELSE mitogenome_data.cox3_trans END,
                cytb_trans = CASE WHEN mitogenome_data.cytb_trans IS NULL THEN EXCLUDED.cytb_trans ELSE mitogenome_data.cytb_trans END,
                nad1_trans = CASE WHEN mitogenome_data.nad1_trans IS NULL THEN EXCLUDED.nad1_trans ELSE mitogenome_data.nad1_trans END,
                nad2_trans = CASE WHEN mitogenome_data.nad2_trans IS NULL THEN EXCLUDED.nad2_trans ELSE mitogenome_data.nad2_trans END,
                nad3_trans = CASE WHEN mitogenome_data.nad3_trans IS NULL THEN EXCLUDED.nad3_trans ELSE mitogenome_data.nad3_trans END,
                nad4_trans = CASE WHEN mitogenome_data.nad4_trans IS NULL THEN EXCLUDED.nad4_trans ELSE mitogenome_data.nad4_trans END,
                nad4l_trans = CASE WHEN mitogenome_data.nad4l_trans IS NULL THEN EXCLUDED.nad4l_trans ELSE mitogenome_data.nad4l_trans END,
                nad5_trans = CASE WHEN mitogenome_data.nad5_trans IS NULL THEN EXCLUDED.nad5_trans ELSE mitogenome_data.nad5_trans END,
                nad6_trans = CASE WHEN mitogenome_data.nad6_trans IS NULL THEN EXCLUDED.nad6_trans ELSE mitogenome_data.nad6_trans END
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
                "seq_date": row_dict["seq_date"],
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

            cursor.execute(upsert_query, params)
            returned = cursor.fetchone()
            conn.commit()

            print("📌 Final stored values:")
            print(dict(zip([
                "annotation", "extra_genes", "missing_genes", "order_correct", "passed", "length_emma", "seqlength_12s",
                "seqlength_16s", "seqlength_co1", "cds_no", "trna_no", "rrna_no", "rrna12s",
                "rrna16s", "atp6", "atp8", "cox1", "cox2", "cox3", "cytb", "nad1", "nad2", "nad3", "nad4", "nad4l",
                "nad5", "nad6", "trna_phe", "trna_val", "trna_leuuag", "trna_leuuaa", "trna_ile", "trna_met",
                "trna_thr", "trna_pro", "trna_lys", "trna_asp", "trna_glu", "trna_sergcu", "trna_seruga",
                "trna_tyr", "trna_cys", "trna_trp", "trna_ala", "trna_asn", "trna_gly", "trna_arg", "trna_his",
                "trna_gln", "atp6_trans", "atp8_trans", "cox1_trans", "cox2_trans", "cox3_trans", "cytb_trans", "nad1_trans",
                "nad2_trans", "nad3_trans", "nad4_trans", "nad4l_trans", "nad5_trans", "nad6_trans"
            ], returned)))

            # Field-by-field comparison
            field_names = [
                "annotation", "extra_genes", "missing_genes", "order_correct", "passed", "length_emma", "seqlength_12s",
                "seqlength_16s", "seqlength_co1", "cds_no", "trna_no", "rrna_no", "rrna12s",
                "rrna16s", "atp6", "atp8", "cox1", "cox2", "cox3", "cytb", "nad1", "nad2", "nad3", "nad4", "nad4l",
                "nad5", "nad6", "trna_phe", "trna_val", "trna_leuuag", "trna_leuuaa", "trna_ile", "trna_met",
                "trna_thr", "trna_pro", "trna_lys", "trna_asp", "trna_glu", "trna_sergcu", "trna_seruga",
                "trna_tyr", "trna_cys", "trna_trp", "trna_ala", "trna_asn", "trna_gly", "trna_arg", "trna_his",
                "trna_gln", "atp6_trans", "atp8_trans", "cox1_trans", "cox2_trans", "cox3_trans", "cytb_trans", "nad1_trans",
                "nad2_trans", "nad3_trans", "nad4_trans", "nad4l_trans", "nad5_trans", "nad6_trans"
            ]
            preserved_fields = []

            for idx, field in enumerate(field_names):
                if params[field] is None:
                    continue  # Skipped intentionally
                elif returned[idx] != params[field]:
                    preserved_fields.append(field)

            if preserved_fields:
                print(f"⚠️ Existing values preserved for: {', '.join(preserved_fields)}")
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
