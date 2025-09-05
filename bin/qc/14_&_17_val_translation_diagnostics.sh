#!/bin/bash

# Summarise .val translation-related errors across GenBank staging folders
STAGING_ROOT="/scratch/pawsey0964/tpeirce/_MITOGENOMES/TREND_remaining/genbank_staging"
OUTFILE="output/val_translation_issues_summary.csv"
PASSFAIL="output/val_passfail_summary.csv"
echo "og_id,seqid,feature,error_type,detail" > "$OUTFILE"
echo "og_id,seqid,status" > "$PASSFAIL"

find "$STAGING_ROOT" -type f -name "*.val" | while read -r valfile; do
  dir=$(dirname "$valfile")
  seqid=$(basename "$(dirname "$dir")")
  og_id=$(basename "$(dirname "$(dirname "$dir")")")

  # Extract error lines and write to full error report
  errors=$(grep -E "InternalStop|MisMatchAA|StartCodon|NoStop" "$valfile")

  if [[ -n "$errors" ]]; then
    echo "$errors" | while read -r line; do
      feature=$(echo "$line" | sed -n 's/.*FEATURE: CDS: \(.*\)\s\{1,\}\[.*/\1/p')
      error_type=$(echo "$line" | sed -n 's/^.*\[SEQ_FEAT\.\(.*\)\].*$/\1/p')
      detail=$(echo "$line" | sed -n 's/^.*\] //p')

      echo "$og_id,$seqid,$feature,$error_type,$detail" >> "$OUTFILE"
    done
    echo "$og_id,$seqid,FAIL" >> "$PASSFAIL"
  else
    echo "$og_id,$seqid,PASS" >> "$PASSFAIL"
  fi

done

echo "✅ CDS translation error summary written to $OUTFILE"
echo "✅ Pass/fail summary written to $PASSFAIL"
