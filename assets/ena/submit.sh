#!/usr/bin/env bash
# Post a Webin drop-box submission and save the receipt.
#
#   ./submit.sh <test|prod> <label> <submission.xml> [project.xml]
#
# Receipts land in receipts/<service>_<label>_<timestamp>.xml. Credentials come from
# the export lines of ~/nextflow_secrets only -- the file's first two lines invoke
# `nextflow secrets set`, so it must never be sourced wholesale.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SECRETS="/home/tpeirce/nextflow_secrets"

service="${1:?service must be test or prod}"
label="${2:?label required, e.g. add|release|modify}"
submission="${3:?submission xml required}"
payload="${4:-}"

case "$service" in
    test) url="https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/" ;;
    prod) url="https://www.ebi.ac.uk/ena/submit/drop-box/submit/" ;;
    *)    echo "service must be 'test' or 'prod', got '$service'" >&2; exit 2 ;;
esac

eval "$(grep -E '^export WEBIN_(USER|PASS)=' "$SECRETS")"
: "${WEBIN_USER:?WEBIN_USER not found in $SECRETS}"
: "${WEBIN_PASS:?WEBIN_PASS not found in $SECRETS}"

mkdir -p "$HERE/receipts"
receipt="$HERE/receipts/${service}_${label}_$(date +%Y%m%dT%H%M%S).xml"

args=(-s --max-time 300 -u "$WEBIN_USER:$WEBIN_PASS" -F "SUBMISSION=@$submission")
# The payload part name is the record type; PROJECT is the only one this workflow submits.
[ -n "$payload" ] && args+=(-F "PROJECT=@$payload")

echo "POST $url  (submission=$submission${payload:+, project=$payload})"
curl "${args[@]}" "$url" > "$receipt"
echo "receipt: $receipt"

python3 - "$receipt" <<'PY'
import sys, xml.etree.ElementTree as ET

path = sys.argv[1]
try:
    root = ET.parse(path).getroot()
except ET.ParseError as exc:
    sys.exit(f"receipt is not valid XML ({exc}); raw body left at {path}")

ok = root.get("success") == "true"
print(f"success={root.get('success')}")

for project in root.findall("PROJECT"):
    print(f"  PROJECT alias={project.get('alias')} accession={project.get('accession')}")
for tag in ("ERROR", "INFO"):
    for msg in root.iter(tag):
        print(f"  {tag}: {(msg.text or '').strip()}")

sys.exit(0 if ok else 1)
PY
