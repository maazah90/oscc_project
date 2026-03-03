#!/bin/bash

BASE=/media/maazah/Expansion/oscc_project
SRR=$1

# Check if SRR is provided
if [ -z "$SRR" ]; then
  echo "Usage: $0 <SRR_ACCESSION>"
  exit 1
fi

echo "==============================="
echo "Processing $SRR"
echo "==============================="

# Retry loop for downloading .sra file with prefetch
echo "Starting prefetch download..."
until prefetch --transport http $SRR -O $BASE/sra; do
  echo "Download failed for $SRR. Retrying in 30 seconds..."
  sleep 30
done

echo "$SRR .sra download successful."

# Retry loop for fasterq-dump conversion
echo "Starting fasterq-dump conversion..."
until fasterq-dump "$SRA_PATH" \
  --split-3 \
  --threads 3 \
  --temp "$BASE/tmp" \
  --outdir "$BASE/fastq" \
  --disable-multithreading; do
  echo "Conversion failed. Retrying in 30 seconds..."
  sleep 30
done

echo "$SRR FASTQ conversion successful."

# Compress FASTQ files
echo "Compressing FASTQ files..."
gzip $BASE/fastq/${SRR}_*.fastq

# Remove SRA file to save space
echo "Cleaning up SRA file..."
rm -rf $BASE/sra/$SRR

echo "$SRR processing complete."#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <SRR_ACCESSION>"
  exit 1
fi

SRR=$1
BASE=/media/maazah/Expansion/oscc_project

echo "Processing $SRR"

prefetch $SRR -O $BASE/sra

fasterq-dump $BASE/sra/$SRR/$SRR.sra \
  --split-3 \
  --threads 3 \
  --temp $BASE/tmp \
  --outdir $BASE/fastq

gzip $BASE/fastq/${SRR}*.fastq

rm -rf $BASE/sra/$SRR

echo "FASTQ complete"
