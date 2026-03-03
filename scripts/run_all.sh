#!/bin/bash

BASE=/media/maazah/Expansion/oscc_project
LIST=srr_list.txt
SCRIPT=./download_and_convert.sh   # <-- change this to the name of your first script

echo "======================================"
echo "Starting batch processing"
echo "======================================"

while read SRR; do

  # Skip empty lines
  [ -z "$SRR" ] && continue

  FASTQ1="$BASE/fastq/${SRR}_1.fastq.gz"
  FASTQ2="$BASE/fastq/${SRR}_2.fastq.gz"

  echo "--------------------------------------"
  echo "Checking $SRR"
  echo "--------------------------------------"

  # Skip if already completed
  if [ -f "$FASTQ1" ] || [ -f "$FASTQ2" ]; then
    echo "$SRR already completed. Skipping."
    continue
  fi

  echo "Processing $SRR ..."
  
  # Run your single-SRR script with logging
  until bash "$SCRIPT" "$SRR" > "$BASE/logs/${SRR}.log" 2>&1; do
    echo "$SRR failed. Retrying in 60 seconds..."
    sleep 60
  done

  echo "$SRR finished successfully."

done < "$LIST"

echo "======================================"
echo "Batch processing complete"
echo "======================================"
