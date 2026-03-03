#!/bin/bash

# Usage:
# bash align_wes.sh SRR31443518

set -e
set -o pipefail

SAMPLE=$1

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample name provided."
    echo "Usage: bash align_wes.sh SAMPLE_ID"
    exit 1
fi

BASE=/media/maazah/Expansion/oscc_project
REF=$BASE/reference/GRCh38.fa
PICARD=$BASE/reference/picard.jar

FASTQ1=$BASE/fastq/${SAMPLE}_1.fastq.gz
FASTQ2=$BASE/fastq/${SAMPLE}_2.fastq.gz

OUTDIR=$BASE/aligned
TMP=$BASE/tmp

mkdir -p $OUTDIR
mkdir -p $TMP

echo "==============================="
echo "Processing $SAMPLE"
echo "==============================="

# 1️⃣ BWA MEM alignment
bwa mem -t 4 \
-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
$REF \
$FASTQ1 \
$FASTQ2 | \
samtools view -b - > $OUTDIR/${SAMPLE}.bam

# 2️⃣ Sort BAM (controlled memory usage)
samtools sort -@ 4 -m 3G \
-o $OUTDIR/${SAMPLE}.sorted.bam \
$OUTDIR/${SAMPLE}.bam

rm $OUTDIR/${SAMPLE}.bam

# 3️⃣ Mark duplicates (Picard)
java -jar $PICARD MarkDuplicates \
I=$OUTDIR/${SAMPLE}.sorted.bam \
O=$OUTDIR/${SAMPLE}.dedup.bam \
M=$OUTDIR/${SAMPLE}.metrics.txt \
CREATE_INDEX=false \
VALIDATION_STRINGENCY=SILENT

rm $OUTDIR/${SAMPLE}.sorted.bam

# 4️⃣ Index final BAM
samtools index $OUTDIR/${SAMPLE}.dedup.bam

echo "==============================="
echo "$SAMPLE alignment complete."
echo "Output: $OUTDIR/${SAMPLE}.dedup.bam"
echo "==============================="
