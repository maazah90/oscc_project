#!/bin/bash

set -e
set -o pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

REF="$PROJECT_DIR/reference/GRCh38.fa"
DBSNP="$PROJECT_DIR/known_sites/dbsnp_146.vcf.gz"
MILLS="$PROJECT_DIR/known_sites/Mills_and_1000G_gold_standard.indels.vcf.gz"
BAMDIR="$PROJECT_DIR/bam"
RESULTS="$PROJECT_DIR/results/recal_trimmed"
LOGS="$PROJECT_DIR/logs_results"
INTERVALS="$PROJECT_DIR/intervals/exome_targets.bed"

mkdir -p "$RESULTS"
mkdir -p "$LOGS"

SAMPLES=($(cat "$PROJECT_DIR/samples.txt"))
TOTAL=${#SAMPLES[@]}
BATCH=2

for ((i=0;i<$TOTAL;i+=BATCH)); do
    echo "===================================="
    echo "Processing batch: ${SAMPLES[@]:i:BATCH}"
    echo "===================================="

    for SAMPLE in "${SAMPLES[@]:i:BATCH}"; do
    (
        RECAL_BAM="$RESULTS/${SAMPLE}_recal.bam"

        if [ -f "$RECAL_BAM" ]; then
            echo "$SAMPLE already recalibrated → skipping"
            exit 0
        fi

        if [ ! -f "$BAMDIR/${SAMPLE}.dedup.bam" ]; then
            echo "Error: BAM not found for $SAMPLE"
            exit 1
        fi

        # Ensure BAM is indexed
        if [ ! -f "$BAMDIR/${SAMPLE}.dedup.bam.bai" ]; then
            samtools index "$BAMDIR/${SAMPLE}.dedup.bam"
        fi

        echo "Running BQSR for $SAMPLE"

        gatk --java-options "-Xmx6G" BaseRecalibrator \
            -R "$REF" \
            -I "$BAMDIR/${SAMPLE}.dedup.bam" \
            --known-sites "$DBSNP" \
            --known-sites "$MILLS" \
            -L "$INTERVALS" \
            -O "$RESULTS/${SAMPLE}_recal.table" \
            2> "$LOGS/${SAMPLE}_bqsr.log"

        echo "Applying BQSR for $SAMPLE"

        gatk --java-options "-Xmx6G" ApplyBQSR \
            -R "$REF" \
            -I "$BAMDIR/${SAMPLE}.dedup.bam" \
            --bqsr-recal-file "$RESULTS/${SAMPLE}_recal.table" \
            -O "$RECAL_BAM" \
            2>> "$LOGS/${SAMPLE}_bqsr.log"

        samtools index "$RECAL_BAM"

        echo "$SAMPLE completed"

    ) &
    done

    wait

    echo ""
    read -p "Batch finished. Press ENTER for next batch..."

done

echo "BQSR step finished for all samples"
