#!/bin/bash

set -euo pipefail

SAMPLE=$1
THREADS=${2:-2}

if [ -z "$SAMPLE" ]; then
    echo "Usage: bash align_wes.sh SAMPLE_ID THREADS"
    exit 1
fi

BASE="/media/maazah/Expansion/oscc_project"
REF="$BASE/reference/GRCh38.fa"
PICARD="$BASE/reference/picard.jar"

FASTQ1="$BASE/trimmed_fastq/${SAMPLE}_1_paired.fastq.gz"
FASTQ2="$BASE/trimmed_fastq/${SAMPLE}_2_paired.fastq.gz"

OUTDIR="$BASE/bam"
LOGDIR="$BASE/logs_trimmed"   # 👈 updated log folder
mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

LOGFILE="$LOGDIR/alignment_runtime.log"

# Check FASTQ files exist
if [[ ! -f "$FASTQ1" || ! -f "$FASTQ2" ]]; then
    echo "FASTQ files missing for $SAMPLE"
    echo "Expected:"
    echo "$FASTQ1"
    echo "$FASTQ2"
    exit 1
fi

# Skip if already completed
if [ -f "$OUTDIR/${SAMPLE}.dedup.bam" ]; then
    echo "$SAMPLE already completed. Skipping."
    exit 0
fi

START_TIME=$(date +%s)
START_HUMAN=$(date)

echo "=================================="
echo "Processing $SAMPLE"
echo "Start: $START_HUMAN"
echo "Using FASTQ1: $FASTQ1"
echo "Using FASTQ2: $FASTQ2"
echo "=================================="

{
# 1️⃣ BWA MEM
bwa mem -t "$THREADS" \
-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
"$REF" \
"$FASTQ1" \
"$FASTQ2" | \
samtools view -@ "$THREADS" -b - > "$OUTDIR/${SAMPLE}.bam"

# 2️⃣ Sort
samtools sort -@ "$THREADS" -m 2G \
-o "$OUTDIR/${SAMPLE}.sorted.bam" \
"$OUTDIR/${SAMPLE}.bam"

rm "$OUTDIR/${SAMPLE}.bam"

# 3️⃣ Mark duplicates
java -Xmx6g -jar "$PICARD" MarkDuplicates \
I="$OUTDIR/${SAMPLE}.sorted.bam" \
O="$OUTDIR/${SAMPLE}.dedup.bam" \
M="$OUTDIR/${SAMPLE}.metrics.txt" \
CREATE_INDEX=false \
VALIDATION_STRINGENCY=SILENT

rm "$OUTDIR/${SAMPLE}.sorted.bam"

# 4️⃣ Index
samtools index "$OUTDIR/${SAMPLE}.dedup.bam"

} && STATUS="SUCCESS" || STATUS="FAILED"

END_TIME=$(date +%s)
END_HUMAN=$(date)
RUNTIME=$((END_TIME - START_TIME))

echo "=================================="
echo "$SAMPLE finished"
echo "End: $END_HUMAN"
echo "Runtime: $RUNTIME seconds"
echo "Status: $STATUS"
echo "=================================="

# Append to master log
echo -e "$SAMPLE\t$START_HUMAN\t$END_HUMAN\t$RUNTIME\t$STATUS" >> "$LOGFILE"
