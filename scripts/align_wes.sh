#!/bin/bash

set -e
set -o pipefail

SAMPLE=$1
THREADS=$2

if [ -z "$SAMPLE" ]; then
    echo "Usage: bash align_wes.sh SAMPLE_ID THREADS"
    exit 1
fi

if [ -z "$THREADS" ]; then
    THREADS=2
fi

BASE=/media/maazah/Expansion/oscc_project
REF=$BASE/reference/GRCh38.fa
PICARD=$BASE/reference/picard.jar

FASTQ1=$BASE/fastq/${SAMPLE}_1.fastq.gz
FASTQ2=$BASE/fastq/${SAMPLE}_2.fastq.gz

OUTDIR=$BASE/aligned
LOGDIR=$BASE/logs
mkdir -p $OUTDIR
mkdir -p $LOGDIR

LOGFILE=$LOGDIR/alignment_runtime.log

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
echo "=================================="

{
# 1️⃣ BWA MEM
bwa mem -t $THREADS \
-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
$REF \
$FASTQ1 \
$FASTQ2 | \
samtools view -b - > $OUTDIR/${SAMPLE}.bam

# 2️⃣ Sort
samtools sort -@ $THREADS -m 2G \
-o $OUTDIR/${SAMPLE}.sorted.bam \
$OUTDIR/${SAMPLE}.bam

rm $OUTDIR/${SAMPLE}.bam

# 3️⃣ Mark duplicates
java -Xmx4g -jar $PICARD MarkDuplicates \
I=$OUTDIR/${SAMPLE}.sorted.bam \
O=$OUTDIR/${SAMPLE}.dedup.bam \
M=$OUTDIR/${SAMPLE}.metrics.txt \
CREATE_INDEX=false \
VALIDATION_STRINGENCY=SILENT

rm $OUTDIR/${SAMPLE}.sorted.bam

# 4️⃣ Index
samtools index $OUTDIR/${SAMPLE}.dedup.bam

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
echo -e "$SAMPLE\t$START_HUMAN\t$END_HUMAN\t$RUNTIME\t$STATUS" >> $LOGFILE
