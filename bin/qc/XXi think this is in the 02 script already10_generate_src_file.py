import pandas as pd
import os
import re
'''
Using: 
    singularity run $SING/psycopg2:0.1.sif python 10_generate_src_file.py

'''
# Input metadata CSV (output of previous cleaning steps)
metadata_file = "output/bankit_metadata_latlon_cleaned.csv"

# Define root directory where GenBank folders are stored
GENBANK_ROOT = "../genbank_staging"  # <-- update this path accordingly

# Load metadata
df = pd.read_csv(metadata_file)

# Explode rows where multiple SeqIDs are listed (comma-separated)
df["SeqID"] = df["SeqID"].astype(str)
df = df.assign(SeqID=df["SeqID"].str.split(",")).explode("SeqID")
df["SeqID"] = df["SeqID"].str.strip()

# Loop over rows and write a .src file per sample
for _, row in df.iterrows():
    full_seqid = row["SeqID"]
    og_id = row["isolate"].strip()
    print(f"full_seqid: {full_seqid}")
    # Take the first 4 dot-separated parts to match directory name
    base_seqid_parts = full_seqid.split(".")
    base_seqid = ".".join(base_seqid_parts[:4])

    genbank_dir = os.path.join(GENBANK_ROOT, og_id, base_seqid, "genbank")
    if not os.path.isdir(genbank_dir):
        print(f"❌ GenBank folder missing: {genbank_dir}")
        continue

    src_path = os.path.join(genbank_dir, f"{full_seqid}.src")

    with open(src_path, "w") as f:
        f.write("SeqID\tisolate\ttissue_type\tcountry\tlat_lon\tCollection_date\n")
        f.write(f"{full_seqid}\t{row['isolate']}\t{row['tissue_type']}\t{row['country']}\t{row['formatted_lat_lon']}\t{row['Collection_date']}\n")

    print(f"✅ Wrote: {src_path}")
