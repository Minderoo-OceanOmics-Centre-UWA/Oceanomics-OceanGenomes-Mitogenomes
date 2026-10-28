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
# Columns whose presence makes an incoming annotation row "complete" enough to
# supersede a stored one. These are fields annotation_stats.py always writes on a
# successful run (it emits 'ok'/'none' rather than a blank for a clean sample,
# precisely so blank keeps meaning "not computed"), so an absent or blank value
# here means the annotation step did not finish, not that it found nothing.
# Bumped whenever annotation_stats.py changes what it computes, so a stored row
# can be identified as the product of older code rather than silently trusted.
ANNOTATION_STATS_VERSION = "2.0.0"

REQUIRED_FOR_COMPLETE = (
    "passed",
    "completeness_profile",
    "annotation_integrity_status",
)


def is_complete_annotation(params):
    """True when the incoming row carries a finished annotation verdict."""
    for col in REQUIRED_FOR_COMPLETE:
        value = params.get(col)
        if value is None or str(value).strip() == "":
            return False
    return True


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
            "annotation", "extra_genes", "missing_genes", "trna_advisory", "order_correct", "order_variant", "order_variant_taxon_check", "order_deviation",
            "annotation_gaps", "passed", "length_emma", "seqlength_12s",
            "seqlength_16s", "seqlength_co1", "cds_no", "trna_no", "rrna_no", "rrna12s",
            "rrna16s", "atp6", "atp8", "cox1", "cox2", "cox3", "cytb", "nad1", "nad2", "nad3", "nad4", "nad4l",
            "nad5", "nad6", "trna_phe", "trna_val", "trna_leuuag", "trna_leuuaa", "trna_ile", "trna_met",
            "trna_thr", "trna_pro", "trna_lys", "trna_asp", "trna_glu", "trna_sergcu", "trna_seruga",
            "trna_tyr", "trna_cys", "trna_trp", "trna_ala", "trna_asn", "trna_gly", "trna_arg", "trna_his",
            "trna_gln", "atp6_trans", "atp8_trans", "cox1_trans", "cox2_trans", "cox3_trans", "cytb_trans", "nad1_trans",
            "nad2_trans", "nad3_trans", "nad4_trans", "nad4l_trans", "nad5_trans", "nad6_trans",
            # Added by migration 026. Every one of these is refreshed on each
            # complete upload: a column present in the INSERT list but absent from
            # the DO UPDATE list is exactly how a value goes stale, and
            # tests/unit/test_push_emma_column_parity.py asserts the two agree.
            "order_status", "order_deviation_detail", "expected_lineage_genes",
            "missing_expected_lineage_genes", "lineage_gene_advisories",
            "mitos_unmapped_features", "mitos_known_auxiliary_features", "duplicate_loci",
            "annotation_integrity_status", "annotation_integrity_issues",
            "atp9", "atp9_trans", "mtmuts", "mtmuts_trans",
            # Added by genomes_db_review migration 068: the gene-set standard `passed`
            # and missing_genes were judged against. REQUIRED_FOR_COMPLETE already
            # demands it; until 068 it was then thrown away.
            "completeness_profile"
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
            -- annotation_updated_at is now(), the DB clock rather than the client's.
            -- "Was this row refreshed by that run" must not depend on a worker node's
            -- time. annotation_stats_version comes from the CSV.
            --
            -- og_num is a generated column maintained by the DB (SUBSTRING(og_id FROM 3)).
            -- It is intentionally absent from this column list: naming it here would error.
            INSERT INTO mitogenome_data (
                og_id, tech, seq_date, code, annotation, extra_genes, missing_genes, trna_advisory, order_correct, order_variant, order_variant_taxon_check, order_deviation, annotation_gaps,
                passed, length_emma, seqlength_12s,
                seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, rrna12s,
                rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, 
                nad5, nad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, 
                trna_thr, trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, 
                trna_tyr, trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, 
                trna_gln, atp6_trans, atp8_trans, cox1_trans, cox2_trans, cox3_trans, cytb_trans, nad1_trans, 
                nad2_trans, nad3_trans, nad4_trans, nad4l_trans, nad5_trans, nad6_trans,
                order_status, order_deviation_detail, expected_lineage_genes,
                missing_expected_lineage_genes, lineage_gene_advisories,
                mitos_unmapped_features, mitos_known_auxiliary_features, duplicate_loci,
                annotation_integrity_status, annotation_integrity_issues,
                atp9, atp9_trans, mtmuts, mtmuts_trans, completeness_profile,
                annotation_stats_version, annotation_updated_at
            )
            VALUES (
                %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s, %(extra_genes)s, %(missing_genes)s, %(trna_advisory)s, %(order_correct)s, %(order_variant)s, %(order_variant_taxon_check)s, %(order_deviation)s,
                %(annotation_gaps)s, %(passed)s, %(length_emma)s, %(seqlength_12s)s,
                %(seqlength_16s)s, %(seqlength_co1)s, %(cds_no)s, %(trna_no)s, %(rrna_no)s, %(rrna12s)s,
                %(rrna16s)s, %(atp6)s, %(atp8)s, %(cox1)s, %(cox2)s, %(cox3)s, %(cytb)s, %(nad1)s, %(nad2)s, %(nad3)s, %(nad4)s, %(nad4l)s,
                %(nad5)s, %(nad6)s, %(trna_phe)s, %(trna_val)s, %(trna_leuuag)s, %(trna_leuuaa)s, %(trna_ile)s, %(trna_met)s,
                %(trna_thr)s, %(trna_pro)s, %(trna_lys)s, %(trna_asp)s, %(trna_glu)s, %(trna_sergcu)s, %(trna_seruga)s,
                %(trna_tyr)s, %(trna_cys)s, %(trna_trp)s, %(trna_ala)s, %(trna_asn)s, %(trna_gly)s, %(trna_arg)s, %(trna_his)s,
                %(trna_gln)s, %(atp6_trans)s, %(atp8_trans)s, %(cox1_trans)s, %(cox2_trans)s, %(cox3_trans)s, %(cytb_trans)s, %(nad1_trans)s, %(nad2_trans)s, 
                %(nad3_trans)s, %(nad4_trans)s, %(nad4l_trans)s, %(nad5_trans)s, %(nad6_trans)s,
                %(order_status)s, %(order_deviation_detail)s, %(expected_lineage_genes)s,
                %(missing_expected_lineage_genes)s, %(lineage_gene_advisories)s,
                %(mitos_unmapped_features)s, %(mitos_known_auxiliary_features)s, %(duplicate_loci)s,
                %(annotation_integrity_status)s, %(annotation_integrity_issues)s,
                %(atp9)s, %(atp9_trans)s, %(mtmuts)s, %(mtmuts_trans)s, %(completeness_profile)s,
                %(annotation_stats_version)s, now()
            )
            ON CONFLICT (og_id, tech, seq_date, code) 
            DO UPDATE SET
                annotation = EXCLUDED.annotation,
                extra_genes = EXCLUDED.extra_genes,
                missing_genes = EXCLUDED.missing_genes,
                trna_advisory = EXCLUDED.trna_advisory,
                order_correct = EXCLUDED.order_correct,
                order_variant = EXCLUDED.order_variant,
                order_variant_taxon_check = EXCLUDED.order_variant_taxon_check,
                order_deviation = EXCLUDED.order_deviation,
                annotation_gaps = EXCLUDED.annotation_gaps,
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
                nad6_trans = EXCLUDED.nad6_trans,
                order_status = EXCLUDED.order_status,
                order_deviation_detail = EXCLUDED.order_deviation_detail,
                expected_lineage_genes = EXCLUDED.expected_lineage_genes,
                missing_expected_lineage_genes = EXCLUDED.missing_expected_lineage_genes,
                lineage_gene_advisories = EXCLUDED.lineage_gene_advisories,
                mitos_unmapped_features = EXCLUDED.mitos_unmapped_features,
                mitos_known_auxiliary_features = EXCLUDED.mitos_known_auxiliary_features,
                duplicate_loci = EXCLUDED.duplicate_loci,
                annotation_integrity_status = EXCLUDED.annotation_integrity_status,
                annotation_integrity_issues = EXCLUDED.annotation_integrity_issues,
                atp9 = EXCLUDED.atp9,
                atp9_trans = EXCLUDED.atp9_trans,
                mtmuts = EXCLUDED.mtmuts,
                mtmuts_trans = EXCLUDED.mtmuts_trans,
                completeness_profile = EXCLUDED.completeness_profile,
                annotation_stats_version = EXCLUDED.annotation_stats_version,
                annotation_updated_at = EXCLUDED.annotation_updated_at
            RETURNING annotation, extra_genes, missing_genes, trna_advisory, order_correct, order_variant, order_variant_taxon_check, order_deviation, annotation_gaps,
                passed, length_emma, seqlength_12s,
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
                "trna_advisory": row_dict.get("trna_advisory"),
                "order_correct": row_dict["order_correct"],
                # .get() rather than [] on purpose: an annotation_stats.csv written
                # before these columns existed must still upsert. This is also the
                # one site of the six that would fail SILENTLY on a missed name --
                # the SQL sites raise on an unknown placeholder, this one would just
                # write nothing.
                "order_variant": row_dict.get("order_variant"),
                "order_variant_taxon_check": row_dict.get("order_variant_taxon_check"),
                "order_deviation": row_dict.get("order_deviation"),
                "annotation_gaps": row_dict.get("annotation_gaps"),
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
                # Migration 026. Optional lookups: an annotation CSV produced
                # before these columns existed simply yields None, which is why
                # older CSVs stay uploadable. annotation_stats.py writes an
                # explicit 'ok'/'none' rather than a blank for a clean sample, so
                # NULL here means "produced by older code", not "clean".
                "order_status": row_dict.get("order_status"),
                "order_deviation_detail": row_dict.get("order_deviation_detail"),
                "expected_lineage_genes": row_dict.get("expected_lineage_genes"),
                "missing_expected_lineage_genes": row_dict.get("missing_expected_lineage_genes"),
                "lineage_gene_advisories": row_dict.get("lineage_gene_advisories"),
                "mitos_unmapped_features": row_dict.get("mitos_unmapped_features"),
                "mitos_known_auxiliary_features": row_dict.get("mitos_known_auxiliary_features"),
                "duplicate_loci": row_dict.get("duplicate_loci"),
                "annotation_integrity_status": row_dict.get("annotation_integrity_status"),
                "annotation_integrity_issues": row_dict.get("annotation_integrity_issues"),
                # mitos_to_emma writes the internal key MTMUTS; the INSDC spelling
                # mtMutS is restored at the submission boundary, not here, so the
                # column name stays lowercase like every other gene column.
                "atp9": row_dict.get("atp9"),
                "atp9_trans": row_dict.get("atp9_trans"),
                "mtmuts": row_dict.get("mtmuts"),
                "mtmuts_trans": row_dict.get("mtmuts_trans"),
                "completeness_profile": row_dict.get("completeness_profile") or None,
                "annotation_stats_version": ANNOTATION_STATS_VERSION,
            }

            # The annotation block is written as ONE unit.
            #
            # This used to refuse any write into a row that already had annotation
            # data, so a rerun that FIXED an annotation left the old values in the
            # database unless someone remembered --force. With the integrity and
            # status columns added that is actively misleading: an
            # annotation_integrity_issues string from a previous run would sit next
            # to an annotation that is now clean, and nothing would say so.
            #
            # A COMPLETE incoming row now supersedes the stored block in full, every
            # annotation column overwritten and the status and issue columns reset to
            # their newly computed values. The protection the old rule was reaching
            # for (never mix old and new) is better served by checking that the
            # incoming row is complete than by refusing to write at all, so an
            # INCOMPLETE row still preserves what is there and still needs --force.
            #
            # Same direction as force_db_overwrite defaulting to true: correcting a
            # record is the normal reason to re-run a sample.
            incoming_complete = is_complete_annotation(params)
            if not force_overwrite and not incoming_complete:
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
                        "the incoming annotation row is INCOMPLETE and the stored row "
                        "already has data; pass --force to overwrite anyway."
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
