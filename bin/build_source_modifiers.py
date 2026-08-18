#!/usr/bin/env python3
# Usage:
#   singularity run $SING/psycopg2:0.1.sif python 02_build_source_modifiers.py \
#     --config config.cfg --og-id OG000123 --seq-tech hifi
#
# config.cfg example:
# [postgres]
# dbname=
# user=
# password=
# host=
# port=

import os
import re
import argparse
import configparser
from pathlib import Path

import pandas as pd
import psycopg2

from geo_loc_name_utils import resolve_geo_loc_name, unmapped_warning

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

# ---------------------------
# CLI
# ---------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Build BankIt source modifiers (per-sample) and generate .src files."
    )
    p.add_argument("--config", required=True, help="Path to INI config file with [postgres].")
    p.add_argument("--og-id", required=True, help="Sample OG_ID / isolate.")
    p.add_argument("--seq-tech", required=True, help="Sequencing tech (e.g., hifi, ilmn, hic).")
    p.add_argument("--assembly-id", required=True, help="Full name of assembly (e.g.OG00.ilmn.220022.getorg1770.emma102).")
    return p.parse_args()

# ---------------------------
# Lat/Lon cleanup helpers
# ---------------------------
def _norm_quotes(s: str) -> str:
    if s is None:
        return s
    return (s
            .replace("′", "'")
            .replace("’", "'")
            .replace("'", "'")
            .replace("″", '"')
            .replace("“", '"')
            .replace("”", '"')
            .replace("o", "°")
            .replace("O", "°")
            .strip())

def parse_coordinate(coord_str):
    """
    Accepts DMS like 33° 52' 31.2" S or decimal with/without sign and optional NSEW.
    Returns 'DD.ddddd X' or 'DD.ddddd' if direction cannot be inferred.
    """
    if coord_str is None:
        return None
    coord_str = _norm_quotes(coord_str)
    direction = ""

    mdir = re.search(r"([NSEW])", coord_str, flags=re.IGNORECASE)
    if mdir:
        direction = mdir.group(1).upper()
        coord_str = re.sub(r"[NSEW]", "", coord_str, flags=re.IGNORECASE).strip()

    mdms = re.match(r"(\d+(?:\.\d+)?)\s*[°:\s]\s*(\d+(?:\.\d+)?)\s*['\s]?\s*(\d+(?:\.\d+)?)?", coord_str)
    if mdms:
        deg = float(mdms.group(1))
        minutes = float(mdms.group(2))
        seconds = float(mdms.group(3)) if mdms.group(3) else 0.0
        dec = deg + minutes/60 + seconds/3600
    else:
        try:
            dec = float(coord_str)
        except ValueError:
            return None
        if not direction:
            if -90 <= dec <= 90:
                direction = "S" if dec < 0 else "N"
            else:
                direction = "W" if dec < 0 else "E"

    dec = abs(dec)
    return f"{dec:.5f} {direction}" if direction else f"{dec:.5f}"

def smart_split_latlon(value):
    if value is None:
        return None, None
    value = _norm_quotes(value)
    if value.lower() == "unknown":
        return None, None

    lat_match = re.search(r'([NS]?\s*[\d°\'"\.\-\s]+[NS])', value, re.IGNORECASE)
    lon_match = re.search(r'([EW]?\s*[\d°\'"\.\-\s]+[EW])', value, re.IGNORECASE)
    if lat_match and lon_match:
        return lat_match.group(1).strip(), lon_match.group(1).strip()

    parts = value.split()
    if len(parts) >= 2:
        return parts[0], parts[1]
    return None, None

def fallback_parse_latlon(value):
    if value is None:
        return "unknown"
    value = _norm_quotes(value)
    if value.lower() == "unknown":
        return "unknown"

    parts = re.findall(r"\d+°\s*\d+(?:\.\d+)?", value)
    if len(parts) >= 2:
        lat_dec = parse_coordinate(parts[0] + " S")
        lon_dec = parse_coordinate(parts[1] + " E")
        if lat_dec and lon_dec:
            return f"{lat_dec} {lon_dec}"

    split = value.split()
    if len(split) >= 2:
        lat_dec = parse_coordinate(split[0] + " S")
        lon_dec = parse_coordinate(split[1] + " E")
        if lat_dec and lon_dec:
            return f"{lat_dec} {lon_dec}"

    return "unknown"

def convert_full_latlon(raw_value):
    lat_raw, lon_raw = smart_split_latlon(raw_value)
    lat_dec = parse_coordinate(lat_raw) if lat_raw else None
    lon_dec = parse_coordinate(lon_raw) if lon_raw else None
    if lat_dec and lon_dec:
        return f"{lat_dec} {lon_dec}"
    return fallback_parse_latlon(raw_value)

# ---------------------------
# Source-modifier validation
# ---------------------------
# INSDC lat_lon: a latitude with a N/S hemisphere then a longitude with an E/W
# hemisphere. Values that parse into something else (two latitudes, a missing
# hemisphere, an out-of-range magnitude) reach table2asn as
# SEQ_DESCR.LatLonFormat, and EMBOSS then demotes them to a stray /note on EMBL
# conversion so ENA never sees a lat_lon at all. Catch them here instead.
LATLON_RE = re.compile(r"^(\d+(?:\.\d+)?) ([NS]) (\d+(?:\.\d+)?) ([EW])$")

def valid_lat_lon(value):
    """True when `value` is a well-formed INSDC lat_lon within coordinate range."""
    match = LATLON_RE.match(str(value).strip())
    if not match:
        return False
    lat, _lat_hem, lon, _lon_hem = match.groups()
    return float(lat) <= 90.0 and float(lon) <= 180.0

# ---------------------------
# DB query
# ---------------------------
def fetch_bankit_metadata(conn, og_id, tech):
    query = """
        WITH sequencing_info AS (
            SELECT
                seq.og_id,
                seq.seq_date,
                CASE
                    WHEN seq.sequencing_id IS NOT NULL
                    THEN SPLIT_PART(seq.sequencing_id, '_', 1)
                    ELSE NULL
                END AS tissue_id
            FROM sequencing seq
        )
        SELECT
            STRING_AGG(
                DISTINCT CONCAT(m.og_id, '.', m.tech, '.', m.seq_date, '.', m.code, '.', m.annotation),
                ', '
            ) AS "SeqID",
            s.og_id AS isolate,
            COALESCE(t.tissue, 'Unknown') AS tissue_type,
            CASE
                -- A NULL country used to collapse straight to 'Unknown', discarding
                -- any locality recorded alongside it. Some rows carry no country but
                -- a locality that names one ('Israel, Elat, Gulf of Aquaba'), so pass
                -- the locality through and let resolve_geo_loc_name promote its
                -- leading token when that token is a controlled value.
                WHEN s.country IS NULL OR s.country = '' THEN
                    CASE
                        WHEN s.location IS NULL OR s.location = '' THEN 'Unknown'
                        ELSE s.location
                    END
                WHEN s.state IS NULL OR s.state = '' THEN
                    CASE
                        WHEN s.location IS NULL OR s.location = '' THEN s.country
                        ELSE CONCAT(s.country, ': ', s.location)
                    END
                ELSE
                    CASE
                        WHEN s.location IS NULL OR s.location = '' THEN CONCAT(s.country, ': ', s.state)
                        ELSE CONCAT(s.country, ': ', s.state, ', ', s.location)
                    END
            END AS country,
            COALESCE(
                NULLIF(CONCAT_WS(' ', s.latitude_collection, s.longitude_collection), ''),
                'unknown'
            ) AS lat_lon,
            s.date_collected
        FROM sample s
        LEFT JOIN mitogenome_data m ON s.og_id = m.og_id
        LEFT JOIN sequencing_info si ON m.og_id = si.og_id AND m.seq_date = si.seq_date
        LEFT JOIN tissue t ON si.tissue_id = t.tissue_id
        WHERE m.og_id = %s AND m.tech = %s
        GROUP BY s.og_id, t.tissue, s.country, s.state, s.location, s.latitude_collection, s.longitude_collection, s.date_collected;
    """
    return pd.read_sql_query(query, conn, params=[og_id, tech])

# ---------------------------
# Main
# ---------------------------
def main():
    args = parse_args()

    db_cfg = load_db_config(args.config)

    Path("output").mkdir(parents=True, exist_ok=True)

    with psycopg2.connect(**db_cfg) as conn:
        print(f"🔍 Fetching metadata for sample: {args.og_id} (tech: {args.seq_tech})")
        df = fetch_bankit_metadata(conn, args.og_id, args.seq_tech)

    if df.empty:
        print(f"❌ No metadata found for sample {args.og_id} with tech {args.seq_tech}")
        pd.DataFrame().to_csv(f"{args.og_id}.bankit_metadata.csv", index=False)
        pd.DataFrame().to_csv(f"{args.og_id}.bankit_metadata_latlon_cleaned.csv", index=False)
        return

    # Format dates
    df["Collection_date"] = pd.to_datetime(df["date_collected"], errors="coerce")
    df["Collection_date"] = df["Collection_date"].dt.strftime("%d-%b-%Y")
    df["Collection_date"] = df["Collection_date"].fillna("Unknown")
    df = df.drop(columns=["date_collected"])

    # Write raw metadata
    out_csv = f"{args.og_id}.bankit_metadata.csv"
    df.to_csv(out_csv, index=False)
    print(f"📁 Metadata written to: {out_csv}")

    # Lat/Lon cleanup + write cleaned
    df["formatted_lat_lon"] = df["lat_lon"].apply(convert_full_latlon)
    cleaned_csv = f"{args.og_id}.bankit_metadata_latlon_cleaned.csv"
    df.to_csv(cleaned_csv, index=False)
    print(f"✅ Updated metadata with cleaned lat/lon values: {cleaned_csv}")

    # Expand SeqIDs and write .src files into the working directory
    df_exp = df.copy()
    df_exp["SeqID"] = df_exp["SeqID"].astype(str).str.split(",")
    df_exp = df_exp.explode("SeqID")
    df_exp["SeqID"] = df_exp["SeqID"].str.strip()

    # Keep only the SeqID for the assembly actually being processed. The SQL
    # STRING_AGG returns every assembly variant for this og_id+tech (e.g. both
    # getorg1770 and getorg1770reseed), but table2asn submits one assembly and
    # must receive exactly one .src; staging several trips a "source qualifier
    # file does not exist" fatal. Match on the assembly prefix passed in (the
    # SeqID is "<assembly-id>.<annotation>", so require a trailing dot so that
    # getorg1770 does not also match getorg1770reseed). Fall back to writing all
    # (the previous behaviour) only if nothing matches, so this can never drop
    # more than it should.
    _aid = (args.assembly_id or "").strip()
    if _aid:
        _match = df_exp["SeqID"].apply(lambda s: s == _aid or s.startswith(_aid + "."))
        if _match.any():
            df_exp = df_exp[_match].copy()
        else:
            print(f"⚠️  No SeqID matched assembly-id '{_aid}'; writing all "
                  f"{len(df_exp)} source file(s) unfiltered.", flush=True)

    # Leave lat_lon empty when we have no usable coordinate. NCBI's validator
    # rejects the literal "unknown" (SEQ_DESCR.LatLonFormat); an empty cell in the
    # .src simply omits the modifier, which is the correct way to signal "no value".
    df_exp["formatted_lat_lon"] = df_exp["lat_lon"].apply(convert_full_latlon)
    df_exp["formatted_lat_lon"] = df_exp["formatted_lat_lon"].fillna("")
    _unknown = df_exp["formatted_lat_lon"].str.strip().str.lower().eq("unknown")
    df_exp.loc[_unknown, "formatted_lat_lon"] = ""

    # A coordinate that survived parsing but is not a well-formed INSDC lat_lon
    # (e.g. "25.51000 N 23.43000 S", two latitudes and no longitude) must not be
    # emitted: table2asn flags it and EMBOSS then demotes it to a note, so the
    # bad value reaches ENA disguised as free text. Drop it and say so loudly.
    _bad_latlon = df_exp["formatted_lat_lon"].ne("") & ~df_exp["formatted_lat_lon"].apply(valid_lat_lon)
    for _seqid, _value in zip(df_exp.loc[_bad_latlon, "SeqID"], df_exp.loc[_bad_latlon, "formatted_lat_lon"]):
        print(f"⚠️  {_seqid}: lat_lon '{_value}' is not a valid INSDC coordinate; "
              f"omitting the modifier. Fix the source record.", flush=True)
    df_exp.loc[_bad_latlon, "formatted_lat_lon"] = ""

    # geo_loc_name (the 'country' modifier) must start with a value from the INSDC
    # controlled list, and it is one of only two mandatory fields in ENA checklist
    # ERC000011. The sample table is rebuilt from a spreadsheet, so typos and
    # long-form names ("Kingdom of Tonga", "Austalia") are translated here rather
    # than corrected upstream, where the next refresh would undo the fix.
    # Where the row records no country, a BioSample this project already registered
    # for the same og_id may state one; failing that the INSDC 'not provided' term
    # is written, because geo_loc_name is mandatory and omitting it fails sample
    # registration outright.
    _resolved = [
        resolve_geo_loc_name(_country, _og_id)
        for _country, _og_id in zip(df_exp["country"], df_exp["isolate"])
    ]
    df_exp["country"] = [value for value, _status in _resolved]
    # An unmapped value is left untouched on purpose: table2asn then raises
    # SEQ_DESCR.BadGeoLocNameCode and the validation gate quarantines this one
    # assembly, which is safer than substituting a placeholder and silently
    # replacing the recorded locality.
    for _seqid, (_value, _status) in zip(df_exp["SeqID"], _resolved):
        if _status == "unmapped":
            print(unmapped_warning(_seqid, _value), flush=True)
        elif _status == "derived":
            # An inference from the locality text rather than a recorded country,
            # so say so: it is the one value here that was not asserted by the DB.
            print(f"ℹ️  {_seqid}: no country recorded; derived geo_loc_name "
                  f"{_value!r} from the locality.", flush=True)
        elif _status == "biosample":
            print(f"ℹ️  {_seqid}: no country recorded; recovered geo_loc_name "
                  f"{_value!r} from the registered NCBI BioSample.", flush=True)

    for _, row in df_exp.iterrows():
        full_seqid = row["SeqID"]
        src_path = os.path.join(".", f"{full_seqid}.src")
        with open(src_path, "w") as f:
            f.write("SeqID\tisolate\ttissue_type\tcountry\tlat_lon\tCollection_date\n")
            f.write(
                f"{full_seqid}\t{row['isolate']}\t{row['tissue_type']}\t"
                f"{row['country']}\t{row['formatted_lat_lon']}\t{row['Collection_date']}\n"
            )
        print(f"✅ Wrote: {src_path}")

if __name__ == "__main__":
    main()
