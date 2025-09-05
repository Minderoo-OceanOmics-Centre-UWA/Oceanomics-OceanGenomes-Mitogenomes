# singularity run $SING/psycopg2:0.1.sif python 99_genbank_submitter.py

from datetime import date
import os
import csv
import smtplib
import psycopg2
from pathlib import Path
from email.message import EmailMessage

# ✏️ CONFIGURATION

# Base folder containing nested .sqn files
STAGED_FOLDER = "/scratch/pawsey0964/tpeirce/_MITOGENOMES/TREND_remaining/genbank_staging/_to_submit/g003"
COVER_LETTER_PATH = os.path.join(STAGED_FOLDER, "email_body.txt")
LOG_PATH = os.path.join(STAGED_FOLDER, "submission_log.csv")

# Embargo date if applicable
EMBARGO_DATE = "2027-01-01"  # e.g. "2025-07-15" or None

# Contact details
EMAIL = "tyler.peirce@uwa.edu.au"
YOUR_NAME = "Tyler Peirce"
INSTITUTION = "Minderoo OceanOmics Centre at UWA, University of Western Australia"

# Email settings
SMTP_SERVER = "smtp.yourdomain.edu"  # e.g., smtp.gmail.com
SMTP_PORT = 587
EMAIL_ADDRESS = "tyler.peirce@uwa.edu.au"
EMAIL_PASSWORD = "your_app_password"
EMAIL_RECIPIENT = "gb-sub@ncbi.nlm.nih.gov"
EMAIL_CC = ["tyler.peirce@uwa.edu.au", "oceangenomes@uwa.edu.au"]
DRY_RUN = False  # Set to False to actually send the email

# Database connection
DB_PARAMS = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '203.101.227.69',
    'port': 5432
}


def parse_filename(filename):
    parts = filename.rstrip(".sqn").split(".")
    if len(parts) < 5:
        raise ValueError(f"Filename {filename} does not match expected pattern.")
    return {
        "og_id": parts[0],
        "tech": parts[1],
        "seq_date": parts[2],
        "code": parts[3],
        "annotation": parts[4],
        "filename": filename
    }


def generate_cover_letter(submission_files, output_path, embargo_date=None):
    today = date.today().isoformat()
    body = [
        "Dear GenBank team,\n",
        f"Please find attached {len(submission_files)} mitochondrial genome submission(s) in .sqn format for review and inclusion in GenBank.\n",
        "These submissions were annotated and validated using table2asn.",
        f"Submission date: {today}."
    ]
    if embargo_date:
        body.append(f"Please hold these sequences until {embargo_date}, or until the associated manuscript is published, whichever comes first.")
    else:
        body.append("These sequences may be released immediately upon acceptance.")

    body += [
        "",
        f"Best regards,\n{YOUR_NAME}\n{INSTITUTION}\n{EMAIL}"
    ]

    with open(output_path, "w") as f:
        f.write("\n".join(body))

    print(f"✅ Cover letter written to {output_path}")


def update_submission_log(log_path, submissions, embargo_date=None):
    today = date.today().isoformat()
    fieldnames = ["og_id", "tech", "seq_date", "annotation", "sqn_filename", "submission_date", "embargo_date"]

    file_exists = os.path.isfile(log_path)
    with open(log_path, "a", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        for sub in submissions:
            writer.writerow({
                "og_id": sub["og_id"],
                "tech": sub["tech"],
                "seq_date": sub["seq_date"],
                "annotation": sub["annotation"],
                "sqn_filename": sub["filename"],
                "submission_date": today,
                "embargo_date": embargo_date or ""
            })
    print(f"📊 Submission log updated: {log_path}")


def update_postgres(submissions, db_params):
    conn = psycopg2.connect(**db_params)
    cur = conn.cursor()
    today = date.today().isoformat()
    count = 0

    for sub in submissions:
        cur.execute("""
            UPDATE mitogenome_data
            SET date_submitted_genbank = %s
            WHERE og_id = %s AND tech = %s AND seq_date = %s AND annotation = %s
        """, (today, sub["og_id"], sub["tech"], sub["seq_date"], sub["annotation"]))
        if cur.rowcount > 0:
            count += 1
        else:
            print(f"⚠️ No match found in DB for: {sub['filename']}")

    conn.commit()
    cur.close()
    conn.close()
    print(f"🗂️ {count} row(s) updated in mitogenome_data")


def main():
    staged_files = list(Path(STAGED_FOLDER).rglob("OG*/OG*/genbank/*.sqn"))
    submissions = []

    for file in staged_files:
        try:
            info = parse_filename(file.name)
            submissions.append(info)
        except ValueError as e:
            print(f"Skipping: {e}")

    if not submissions:
        print("❌ No valid .sqn files found.")
        return

    # Create working directory to collect everything
    order_name = Path(STAGED_FOLDER).name
    export_dir = Path(STAGED_FOLDER) / f"{order_name}_submission_package"
    export_dir.mkdir(parents=True, exist_ok=True)

    # Copy .sqn files into working directory
    for file in staged_files:
        dest = export_dir / file.name
        dest.write_bytes(file.read_bytes())

    cover_letter_path = export_dir / "email_body.txt"
    generate_cover_letter(submissions, cover_letter_path, embargo_date=EMBARGO_DATE)
    update_submission_log(LOG_PATH, submissions, embargo_date=EMBARGO_DATE)
    update_postgres(submissions, DB_PARAMS)

    attachments = [str(f) for f in export_dir.glob("*.sqn")]

    instructions_path = export_dir / "email_instructions.txt"
    with open(instructions_path, "w", encoding="utf-8") as f:
        f.write("📬 Manual Email Instructions:\n")
        f.write(f"To: {EMAIL_RECIPIENT}\n")
        f.write(f"CC: {'; '.join(EMAIL_CC)}\n")
        f.write(f"Subject: GenBank mitogenome submission batch {order_name}\n")
        f.write(f"Body: See contents of {cover_letter_path}\n")
        f.write("Attach the following .sqn files:\n")
        for file in attachments:
            f.write(f" - {file}\n")
        f.write("\nThen send the email from your preferred email client.\n")


    print("📬 Manual Email Instructions:")
    print(f"To: {EMAIL_RECIPIENT}")
    print(f"CC: {', '.join(EMAIL_CC)}")
    print("Subject: GenBank mitogenome submission batch")
    print(f"Body: See contents of {cover_letter_path}")
    print("Attach the following .sqn files:")
    for file in attachments:
        print(f" - {file}")
    print("Then send the email from your preferred email client.")
    print(f"📁 Email instructions saved to: {instructions_path}")


if __name__ == "__main__":
    main()
