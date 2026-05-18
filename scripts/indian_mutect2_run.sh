#!/bin/bash

set -e
set -o pipefail

################################
# PROJECT PATH SETUP
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

REF="$PROJECT_DIR/reference/GRCh38.fa"
GNOMAD="$PROJECT_DIR/resources/af-only-gnomad.fixed.vcf.gz"
GNOMAD_COMMON="$PROJECT_DIR/resources/gnomad_common.fixed.vcf.gz"
PON="$PROJECT_DIR/resources/1000g_pon.hg38.fixed.vcf.gz"
INTERVALS="$PROJECT_DIR/intervals/exome_targets.bed"

RESULTS="$PROJECT_DIR/results/mutect2_indian_samples"
LOGS="$PROJECT_DIR/logs_indian_mutect2"

mkdir -p "$RESULTS"
mkdir -p "$RESULTS/tmp"
mkdir -p "$LOGS"

################################
# SAMPLE LIST & BATCH SETUP
################################

SAMPLES=($(cat "$PROJECT_DIR/scripts/indian_srrs_clean.txt"))
TOTAL=${#SAMPLES[@]}
BATCH=2   # 2 samples in parallel

################################
# LOOP OVER BATCHES
################################

for ((i=0;i<$TOTAL;i+=BATCH)); do

    echo "===================================="
    echo "Running batch: ${SAMPLES[@]:i:BATCH}"
    echo "===================================="

    for SAMPLE in "${SAMPLES[@]:i:BATCH}"; do
    (
        echo "------------------------------------"
        echo "Processing $SAMPLE"
        echo "------------------------------------"

        FINAL_VCF="$RESULTS/${SAMPLE}_filtered.vcf.gz"
        UNFILTERED="$RESULTS/${SAMPLE}_unfiltered.vcf.gz"
        RECAL_BAM="$PROJECT_DIR/results/recal_indian/${SAMPLE}_recal.bam"

        if [ -f "$FINAL_VCF" ]; then
            echo "$SAMPLE already complete → skipping"
            exit 0
        fi

        if [ ! -f "$RECAL_BAM" ]; then
            echo "ERROR: $RECAL_BAM not found"
            exit 1
        fi

        if [ ! -f "${RECAL_BAM}.bai" ]; then
            echo "Indexing BAM for $SAMPLE"
            samtools index "$RECAL_BAM"
        fi

        ################################
        # TMP DIR
        ################################
        TMPDIR_SAMPLE="$RESULTS/tmp/$SAMPLE"
        mkdir -p "$TMPDIR_SAMPLE"

        ################################
        # INTERMEDIATE FILES
        ################################
        F1R2="$RESULTS/${SAMPLE}_f1r2.tar.gz"
        PRIORS="$RESULTS/${SAMPLE}_artifact-priors.tar.gz"
        PILEUP="$RESULTS/${SAMPLE}_pileups.table"
        CONTAM="$RESULTS/${SAMPLE}_contamination.table"

        ################################
        # MUTECT2
        ################################

        echo "Running Mutect2 for $SAMPLE"

        gatk --java-options "-Xmx6G -Djava.io.tmpdir=$TMPDIR_SAMPLE" Mutect2 \
            -R "$REF" \
            -I "$RECAL_BAM" \
            -L "$INTERVALS" \
            -tumor "$SAMPLE" \
            --germline-resource "$GNOMAD" \
            --panel-of-normals "$PON" \
            --native-pair-hmm-threads 2 \
            --f1r2-tar-gz "$F1R2" \
            -O "$UNFILTERED" \
            2> "$LOGS/${SAMPLE}_mutect2.log"

        ################################
        # VALIDATION
        ################################

        if [ ! -f "$UNFILTERED" ] || [ ! -f "${UNFILTERED}.stats" ]; then
            echo "ERROR: Mutect2 output missing for $SAMPLE"
            exit 1
        fi

        ################################
        # ORIENTATION BIAS
        ################################

        echo "Learning orientation bias for $SAMPLE"

        gatk --java-options "-Xmx6G -Djava.io.tmpdir=$TMPDIR_SAMPLE" LearnReadOrientationModel \
            -I "$F1R2" \
            -O "$PRIORS"

        ################################
        # CONTAMINATION
        ################################

        echo "Estimating contamination for $SAMPLE"

        gatk --java-options "-Xmx6G -Djava.io.tmpdir=$TMPDIR_SAMPLE" GetPileupSummaries \
            -I "$RECAL_BAM" \
            -V "$GNOMAD_COMMON" \
            -L "$INTERVALS" \
            -O "$PILEUP"

        gatk --java-options "-Xmx6G -Djava.io.tmpdir=$TMPDIR_SAMPLE" CalculateContamination \
            -I "$PILEUP" \
            -O "$CONTAM"

        ################################
        # FILTERING
        ################################

        echo "Filtering variants for $SAMPLE"

        gatk --java-options "-Xmx6G -Djava.io.tmpdir=$TMPDIR_SAMPLE" FilterMutectCalls \
            -R "$REF" \
            -V "$UNFILTERED" \
            --stats "${UNFILTERED}.stats" \
            --contamination-table "$CONTAM" \
            --ob-priors "$PRIORS" \
            -O "$FINAL_VCF"

        ################################
        # INDEX FINAL VCF
        ################################

        gatk IndexFeatureFile -I "$FINAL_VCF"

        ################################
        # CLEANUP
        ################################

        rm -f "$F1R2" "$PRIORS" "$PILEUP" "$CONTAM"

        echo "$SAMPLE finished successfully"

    ) &
    done

    wait

    echo "Batch finished. Continuing..."

done

echo "Pipeline finished"
