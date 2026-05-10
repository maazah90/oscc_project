#!/bin/bash

set -e
set -o pipefail

################################
# PATHS
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

REF="$PROJECT_DIR/reference/GRCh38.fa"
GNOMAD="$PROJECT_DIR/resources/af-only-gnomad.fixed.vcf.gz"
GNOMAD_COMMON="$PROJECT_DIR/resources/gnomad_common.fixed.vcf.gz"

INTERVALS="$PROJECT_DIR/intervals/exome_targets.bed"
BAM_DIR="$PROJECT_DIR/results/recal_trimmed"

RESULTS="$PROJECT_DIR/results/mutect2_paired"
LOGS="$PROJECT_DIR/logs_results"

mkdir -p "$RESULTS" "$RESULTS/tmp" "$LOGS"

################################
# OPTIMISATION SETTINGS
################################

THREADS=4
BATCH=1

################################
# PROCESS FUNCTION
################################

process_pair() {

    TUMOR=$1
    NORMAL=$2

    echo "======================================"
    echo "Processing: $TUMOR vs $NORMAL"
    echo "======================================"

    TUMOR_BAM="$BAM_DIR/${TUMOR}_recal.bam"
    NORMAL_BAM="$BAM_DIR/${NORMAL}_recal.bam"

    FINAL_VCF="$RESULTS/${TUMOR}_vs_${NORMAL}.filtered.vcf.gz"
    UNFILTERED_VCF="$RESULTS/${TUMOR}_vs_${NORMAL}.unfiltered.vcf.gz"

    F1R2="$RESULTS/${TUMOR}_f1r2.tar.gz"
    ORIENT_MODEL="$RESULTS/${TUMOR}_orientation.tar.gz"

    PILEUP="$RESULTS/${TUMOR}_pileup.table"
    CONTAM="$RESULTS/${TUMOR}_contamination.table"

    TMPDIR="$RESULTS/tmp/${TUMOR}"
    mkdir -p "$TMPDIR"

    if [ -f "$FINAL_VCF" ]; then
        echo "Skipping $TUMOR (done)"
        return
    fi

    ################################
    # MUTECT2 (OPTIMISED)
    ################################

    gatk --java-options "-Xmx5G -Djava.io.tmpdir=$TMPDIR" Mutect2 \
        -R "$REF" \
        -I "$TUMOR_BAM" \
        -I "$NORMAL_BAM" \
        -tumor "$TUMOR" \
        -normal "$NORMAL" \
        -L "$INTERVALS" \
        --germline-resource "$GNOMAD" \
        --native-pair-hmm-threads "$THREADS" \
        --f1r2-tar-gz "$F1R2" \
        -O "$UNFILTERED_VCF" \
        2> "$LOGS/${TUMOR}_mutect2.log"

    ################################
    # ORIENTATION MODEL
    ################################

    gatk --java-options "-Xmx2G" LearnReadOrientationModel \
        -I "$F1R2" \
        -O "$ORIENT_MODEL"

    ################################
    # CONTAMINATION
    ################################

    gatk --java-options "-Xmx2G" GetPileupSummaries \
        -I "$TUMOR_BAM" \
        -V "$GNOMAD_COMMON" \
        -L "$INTERVALS" \
        -O "$PILEUP"

    gatk --java-options "-Xmx2G" CalculateContamination \
        -I "$PILEUP" \
        -O "$CONTAM"

    ################################
    # FILTERING
    ################################

    gatk --java-options "-Xmx5G" FilterMutectCalls \
        -R "$REF" \
        -V "$UNFILTERED_VCF" \
        --contamination-table "$CONTAM" \
        --ob-priors "$ORIENT_MODEL" \
        -O "$FINAL_VCF"

    gatk IndexFeatureFile -I "$FINAL_VCF"

    rm -f "$F1R2" "$PILEUP" "$CONTAM"

    echo "Done: $TUMOR vs $NORMAL"
}

################################
# RUN PAIRS SEQUENTIALLY (OPTIMAL MODE)
################################

while read TUMOR NORMAL; do
    process_pair "$TUMOR" "$NORMAL"
done < pairs_srr.txt

echo "All pairs completed"
