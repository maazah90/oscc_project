#!/bin/bash

################################
# USER PATHS
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

REF=$PROJECT_DIR/reference/GRCh38.fa
DBSNP=$PROJECT_DIR/known_sites/dbsnp_146.vcf.gz
MILLS=$PROJECT_DIR/known_sites/Mills_and_1000G_gold_standard.indels.vcf.gz
BAMDIR=$PROJECT_DIR/aligned
RESULTS=$PROJECT_DIR/results
LOGS=$PROJECT_DIR/logs
INTERVALS=$PROJECT_DIR/intervals/exome_targets.bed

# Make sure all result directories exist
mkdir -p "$RESULTS/recal"
mkdir -p "$LOGS"

################################
# SAMPLE LIST & BATCH SETUP
################################
SAMPLES=($(cat $PROJECT_DIR/samples.txt))
TOTAL=${#SAMPLES[@]}
BATCH=2   # safe for 4 cores / 20GB RAM

################################
# LOOP OVER BATCHES
################################
for ((i=0;i<$TOTAL;i+=BATCH)); do
    echo "===================================="
    echo "Processing batch: ${SAMPLES[@]:i:BATCH}"
    echo "===================================="

    for SAMPLE in "${SAMPLES[@]:i:BATCH}"; do

        RECAL_BAM=$RESULTS/recal/${SAMPLE}_recal.bam

        # Skip finished samples
        if [ -f "$RECAL_BAM" ]; then
            echo "$SAMPLE already recalibrated → skipping"
            continue
        fi

        # Check BAM exists
        if [ ! -f "$BAMDIR/${SAMPLE}.dedup.bam" ]; then
            echo "Error: $BAMDIR/${SAMPLE}.dedup.bam not found, skipping $SAMPLE"
            continue
        fi

        echo "Running BQSR for $SAMPLE"

        ################################
        # Base Quality Score Recalibration
        ################################
        gatk --java-options "-Xmx4G" BaseRecalibrator \
            -R $REF \
            -I $BAMDIR/${SAMPLE}.dedup.bam \
            --known-sites $DBSNP \
            --known-sites $MILLS \
            -L $INTERVALS \
            -O $RESULTS/recal/${SAMPLE}_recal.table \
            2> $LOGS/${SAMPLE}_bqsr.log

        if [ $? -ne 0 ]; then
            echo "BaseRecalibrator failed for $SAMPLE"
            continue
        fi

        ################################
        # Apply BQSR
        ################################
        gatk --java-options "-Xmx4G" ApplyBQSR \
            -R $REF \
            -I $BAMDIR/${SAMPLE}.dedup.bam \
            --bqsr-recal-file $RESULTS/recal/${SAMPLE}_recal.table \
            -O $RECAL_BAM \
            2>> $LOGS/${SAMPLE}_bqsr.log

        if [ $? -ne 0 ]; then
            echo "ApplyBQSR failed for $SAMPLE"
            continue
        fi

        ################################
        # Index recalibrated BAM
        ################################
        samtools index $RECAL_BAM
        echo "$SAMPLE BAM recalibrated and indexed"

    done

    echo ""
    read -p "Batch finished. Press ENTER for next batch..."

done

echo "BQSR step finished for all samples"
