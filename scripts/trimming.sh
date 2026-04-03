#!/bin/bash
set -euo pipefail

THREADS=4
TRIMMO="/media/maazah/Expansion/oscc_project/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/media/maazah/Expansion/oscc_project/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

RAW_DIR="/media/maazah/Expansion/oscc_project/fastq"
TRIM_DIR="/media/maazah/Expansion/oscc_project/trimmed_fastq"

mkdir -p "$TRIM_DIR"

while read -r sample; do
  echo "[$(date)] Processing $sample..."

  if [[ ! -f "$RAW_DIR/${sample}_1.fastq.gz" || ! -f "$RAW_DIR/${sample}_2.fastq.gz" ]]; then
    echo "Missing FASTQ files for $sample, skipping..."
    continue
  fi

  java -jar "$TRIMMO" PE -threads "$THREADS" \
    "$RAW_DIR/${sample}_1.fastq.gz" \
    "$RAW_DIR/${sample}_2.fastq.gz" \
    "$TRIM_DIR/${sample}_1_paired.fastq.gz" \
    "$TRIM_DIR/${sample}_1_unpaired.fastq.gz" \
    "$TRIM_DIR/${sample}_2_paired.fastq.gz" \
    "$TRIM_DIR/${sample}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    > "$TRIM_DIR/${sample}_trimmomatic.log" 2>&1

done < samples.txt
