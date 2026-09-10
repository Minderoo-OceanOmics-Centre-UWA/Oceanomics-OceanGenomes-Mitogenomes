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
from datetime import date, datetime
from pathlib import Path

import pandas as pd
import psycopg2

from geo_loc_name_utils import (
    hemispheres_for_geo_loc_name,
    resolve_geo_loc_name,
    unmapped_warning,
)

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

# A hemisphere is information, and where the sample table does not record one it
# must not be invented. Latitudes arrive as unsigned magnitudes -- a Ningaloo
# sample at 22.03 S is stored as "22.03 113.891" -- and defaulting an unsigned
# value to the northern hemisphere put fifteen Western Australian assemblies on
# the wrong side of the equator. table2asn accepts the shape and then rejects the
# value as SEQ_DESCR.LatLonValue because the coordinate contradicts the country,
# which quarantines the whole assembly. So the hemisphere is only ever *read*
# here; where the value carries none it is derived from the resolved
# geo_loc_name, and where the country cannot settle it the modifier is omitted.
# See format_lat_lon below.
AXIS_HEMISPHERES = {"lat": ("N", "S"), "lon": ("E", "W")}

def parse_coordinate(coord_str, axis):
    """Parse one coordinate into (magnitude, hemisphere-or-None).

    `axis` is "lat" or "lon" and fixes which hemisphere letters are legal, so a
    longitude written "23.43 S" is rejected here rather than emitted as a second
    latitude. Accepts DMS like 33° 52' 31.2" S and decimal with or without a sign
    and with or without an NSEW letter.

    Returns None when the text is not a coordinate for this axis at all. Returns
    a None hemisphere when the value simply does not state one, which is a
    different thing and is resolved against the country by the caller.
    """
    if coord_str is None:
        return None
    positive, negative = AXIS_HEMISPHERES[axis]
    coord_str = _norm_quotes(coord_str)
    hemisphere = None

    mdir = re.search(r"([NSEW])", coord_str, flags=re.IGNORECASE)
    if mdir:
        hemisphere = mdir.group(1).upper()
        if hemisphere not in (positive, negative):
            # A latitude labelled E, or a longitude labelled S. No correction is
            # safe to guess, so this is not a coordinate.
            return None
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
        if hemisphere is None and coord_str.strip().startswith(("-", "+")):
            # An explicit sign IS a recorded hemisphere. Only a bare magnitude is
            # silent, and only a bare magnitude falls through to the country.
            hemisphere = negative if dec < 0 else positive

    return abs(dec), hemisphere

def split_lat_lon(value):
    """Split a raw sample-table coordinate into its latitude and longitude text."""
    if value is None:
        return None, None
    value = _norm_quotes(value)
    if value == "" or value.lower() == "unknown":
        return None, None

    lat_match = re.search(r'([NS]?\s*[\d°\'"\.\-\s]+[NS])', value, re.IGNORECASE)
    lon_match = re.search(r'([EW]?\s*[\d°\'"\.\-\s]+[EW])', value, re.IGNORECASE)
    if lat_match and lon_match:
        return lat_match.group(1).strip(), lon_match.group(1).strip()

    # A DMS pair carrying no hemisphere letters ("12° 34.5 113° 20.1"), which the
    # whitespace split below would tear in half.
    dms_parts = re.findall(r"\d+°\s*\d+(?:\.\d+)?", value)
    if len(dms_parts) >= 2:
        return dms_parts[0], dms_parts[1]

    parts = value.split()
    if len(parts) >= 2:
        return parts[0], parts[1]
    return None, None

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

def format_lat_lon(raw_value, geo_loc_name):
    """Render a raw sample-table coordinate as an INSDC lat_lon.

    Returns (value, status). `value` is "" whenever nothing submittable can be
    built, because an empty cell in the .src omits the modifier, and omission is
    the only correct way to say "no value" here: the literal "unknown" is
    rejected as SEQ_DESCR.LatLonFormat, and a guessed hemisphere is rejected as
    SEQ_DESCR.LatLonValue. `status` is one of:

        'ok'            -- a complete, well-formed coordinate
        'absent'        -- nothing recorded; the common case, reported silently
        'unparsable'    -- text that is not a coordinate, or one out of range
        'no_hemisphere' -- a bare magnitude the geo_loc_name cannot place,
                           because the country straddles that axis or is not a
                           mapped country at all
        'conflict'      -- a stated hemisphere the geo_loc_name contradicts,
                           which is precisely what table2asn rejects
    """
    lat_raw, lon_raw = split_lat_lon(raw_value)
    if lat_raw is None and lon_raw is None:
        return "", "absent"

    lat = parse_coordinate(lat_raw, "lat")
    lon = parse_coordinate(lon_raw, "lon")
    if lat is None or lon is None:
        return "", "unparsable"
    lat_dec, lat_hem = lat
    lon_dec, lon_hem = lon

    country_lat_hem, country_lon_hem = hemispheres_for_geo_loc_name(geo_loc_name)

    # A recorded hemisphere the country disagrees with is a bad record, not a
    # value to correct: which of the two is wrong cannot be told from here.
    for stated, implied in ((lat_hem, country_lat_hem), (lon_hem, country_lon_hem)):
        if stated and implied and stated != implied:
            return "", "conflict"

    lat_hem = lat_hem or country_lat_hem
    lon_hem = lon_hem or country_lon_hem
    if not lat_hem or not lon_hem:
        return "", "no_hemisphere"

    formatted = f"{lat_dec:.5f} {lat_hem} {lon_dec:.5f} {lon_hem}"
    return (formatted, "ok") if valid_lat_lon(formatted) else ("", "unparsable")

# Operator-facing reason for each format_lat_lon status that drops a coordinate.
# 'absent' is deliberately not here: most samples record no coordinate at all and
# warning on every one of them would bury the rows that are actually broken.
LATLON_DROP_REASONS = {
    "unparsable": "is not a usable INSDC coordinate",
    "no_hemisphere": ("records no hemisphere, and the geo_loc_name does not imply "
                      "one (add the country to COUNTRY_HEMISPHERE in "
                      "bin/geo_loc_name_utils.py, or sign the latitude at source)"),
    "conflict": "states a hemisphere that its geo_loc_name contradicts",
}

# INSDC collection_date, in the three forms both gates accept: DD-Mmm-YYYY,
# Mmm-YYYY and YYYY. Anything else is SEQ_DESCR.BadCollectionDate at table2asn,
# which the validation gate counts as an ERROR and quarantines the assembly for.
#
# There is deliberately no missing-value term here. 'missing', 'not collected'
# and 'not provided' all pass table2asn, so they look like the obvious fix, but
# they belong to the ENA *sample checklist* vocabulary (ERC000011), not to the
# flatfile qualifier: ENA's own CollectionDateQualifierCheck (sequencetools, the
# validator ena-webin-cli runs) rejects every one of them. Adopting one would
# just move the failure from the first gate to the last. collection_date is not
# a required source qualifier in the ENA flatfile, so where no date is recorded
# the cell is left empty and the modifier is simply omitted -- the same
# resolution valid_lat_lon uses, and the only one that clears both gates.
COLLECTION_DATE_RE = re.compile(
    r"^(?:(?:\d{1,2}-)?(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-)?\d{4}$"
)

def valid_collection_date(value):
    """True when `value` is a well-formed INSDC collection_date that is not in the future."""
    text = str(value).strip()
    if not COLLECTION_DATE_RE.match(text):
        return False
    # Both validators also reject a date later than today (table2asn:
    # "Collection_date is in the future"; ENA: FutureDateException). The sample
    # table carries day/month-transposed rows such as 2026-12-06 that parse
    # cleanly but land in the future, so the format check alone is not enough.
    for fmt in ("%d-%b-%Y", "%b-%Y", "%Y"):
        try:
            parsed = datetime.strptime(text, fmt)
        except ValueError:
            continue
        return parsed.date() <= date.today()
    return False

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
        # Fail rather than write empty tables. The query SELECTs only from sample,
        # but it INNER-filters on mitogenome_data (WHERE m.og_id / m.tech), so an
        # empty result almost always means the mitogenome_data parent row does not
        # exist yet -- not that the specimen has no collection metadata.
        #
        # This used to write two empty CSVs and exit 0, which produced a package
        # with no source modifiers at all and said so only in the task log. That
        # was survivable while submission prep was gated behind the committed
        # upload receipt and the row was therefore guaranteed. Prep now runs with
        # --skip_upload_results, where a sample new to the database genuinely has
        # no row, and a silent empty .src there would reach a submitter looking
        # exactly like a specimen with nothing recorded about it.
        raise SystemExit(
            f"❌ No metadata found for {args.og_id} (tech: {args.seq_tech}). "
            "The mitogenome_data row this query filters on is missing, or the "
            "specimen has no sample record. Upload the assembly results first "
            "(run without --skip_upload_results), or fix the sample record."
        )

    # Format dates
    df["Collection_date"] = pd.to_datetime(df["date_collected"], errors="coerce")
    df["Collection_date"] = df["Collection_date"].dt.strftime("%d-%b-%Y")
    # An unrecorded date leaves the cell empty, which omits the modifier. It used
    # to be filled with the literal "Unknown", which table2asn rejects outright as
    # SEQ_DESCR.BadCollectionDate -- see valid_collection_date above for why no
    # missing-value term is substituted instead.
    df["Collection_date"] = df["Collection_date"].fillna("")
    df = df.drop(columns=["date_collected"])

    # Write raw metadata
    out_csv = f"{args.og_id}.bankit_metadata.csv"
    df.to_csv(out_csv, index=False)
    print(f"📁 Metadata written to: {out_csv}")

    # geo_loc_name (the 'country' modifier) must start with a value from the INSDC
    # controlled list, and it is one of only two mandatory fields in ENA checklist
    # ERC000011. The sample table is rebuilt from a spreadsheet, so typos and
    # long-form names ("Kingdom of Tonga", "Austalia") are translated here rather
    # than corrected upstream, where the next refresh would undo the fix.
    # Where the row records no country, a BioSample this project already registered
    # for the same og_id may state one; failing that the INSDC 'not provided' term
    # is written, because geo_loc_name is mandatory and omitting it fails sample
    # registration outright.
    #
    # Resolved before the coordinate below, and not after it as it once was,
    # because a latitude the sample table records without a hemisphere can only
    # be placed against the country -- and only against the canonical spelling of
    # it, so that the alias table is the single place country names are corrected.
    _resolved = [
        resolve_geo_loc_name(_country, _og_id)
        for _country, _og_id in zip(df["country"], df["isolate"])
    ]
    df["country"] = [value for value, _status in _resolved]
    df["_geo_loc_status"] = [status for _value, status in _resolved]

    # Lat/Lon cleanup + write cleaned
    _formatted = [
        format_lat_lon(_raw, _country)
        for _raw, _country in zip(df["lat_lon"], df["country"])
    ]
    df["formatted_lat_lon"] = [value for value, _status in _formatted]
    df["_lat_lon_status"] = [status for _value, status in _formatted]
    cleaned_csv = f"{args.og_id}.bankit_metadata_latlon_cleaned.csv"
    df.drop(columns=["_geo_loc_status", "_lat_lon_status"]).to_csv(cleaned_csv, index=False)
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

    # A coordinate that cannot be turned into a submittable INSDC lat_lon is
    # dropped rather than guessed at, and the reason is named so the source record
    # can be fixed. Emitting a bad one is worse than emitting none: table2asn
    # flags it and EMBOSS then demotes it to a note, so the wrong value reaches
    # ENA disguised as free text. Warned per SeqID here, after the assembly filter
    # above, so the message names the assembly actually being submitted.
    for _seqid, _raw, _status in zip(
        df_exp["SeqID"], df_exp["lat_lon"], df_exp["_lat_lon_status"]
    ):
        if _status in LATLON_DROP_REASONS:
            print(f"⚠️  {_seqid}: lat_lon '{_raw}' {LATLON_DROP_REASONS[_status]}; "
                  f"omitting the modifier. Fix the source record.", flush=True)

    # Same treatment for a date that is present but not submittable -- a future
    # date, or anything that is not one of the three INSDC forms. Unlike
    # geo_loc_name this is not mandatory anywhere on the ENA flatfile route, so an
    # omitted date costs nothing, whereas a bad one fails the gate for the whole
    # assembly. Most samples reaching here simply have no date recorded at all;
    # that is already an empty cell and is passed over silently rather than
    # warned about on every record.
    df_exp["Collection_date"] = df_exp["Collection_date"].fillna("").astype(str).str.strip()
    _bad_date = df_exp["Collection_date"].ne("") & ~df_exp["Collection_date"].apply(valid_collection_date)
    for _seqid, _value in zip(df_exp.loc[_bad_date, "SeqID"], df_exp.loc[_bad_date, "Collection_date"]):
        print(f"⚠️  {_seqid}: Collection_date '{_value}' is not a submittable INSDC "
              f"date; omitting the modifier. Fix the source record.", flush=True)
    df_exp.loc[_bad_date, "Collection_date"] = ""

    # An unmapped country is left untouched on purpose: table2asn then raises
    # SEQ_DESCR.BadGeoLocNameCode and the validation gate quarantines this one
    # assembly, which is safer than substituting a placeholder and silently
    # replacing the recorded locality. Resolution itself happened before the
    # coordinate block above; only the reporting waits until here, so each message
    # names the assembly being submitted rather than the aggregated row.
    for _seqid, _value, _status in zip(
        df_exp["SeqID"], df_exp["country"], df_exp["_geo_loc_status"]
    ):
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
