#!/bin/bash

################################
# PROJECT PATH SETUP
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

REF=$PROJECT_DIR/reference/GRCh38.fa
GNOMAD=$PROJECT_DIR/resources/af-only-gnomad.vcf.gz
GNOMAD_COMMON=$PROJECT_DIR/resources/gnomad_common.vcf.gz
INTERVALS=$PROJECT_DIR/intervals/exome_targets.bed

RESULTS=$PROJECT_DIR/results
LOGS=$PROJECT_DIR/logs

# Create results/logs folders if they don't exist
mkdir -p $RESULTS/mutect2
mkdir -p $LOGS
mkdir -p $RESULTS/tmp

################################
# FORCE ALL TEMP FILES TO EXTERNAL DRIVE
################################

export TMPDIR=$RESULTS/tmp
export _JAVA_OPTIONS="-Djava.io.tmpdir=$TMPDIR"

echo "GATK temporary folder set to: $TMPDIR"

################################
# SAMPLE LIST & BATCH SETUP
################################

SAMPLES=($(cat $PROJECT_DIR/samples.txt))
TOTAL=${#SAMPLES[@]}
BATCH=2   # adjusted for 4 cores / 20GB RAM

################################
# MUTECT2 LOOP (Tumor-only)
################################

for ((i=0;i<$TOTAL;i+=BATCH)); do

    echo "===================================="
    echo "Running batch ${SAMPLES[@]:i:BATCH}"
    echo "===================================="

    for SAMPLE in "${SAMPLES[@]:i:BATCH}"; do

        FINAL_VCF=$RESULTS/mutect2/${SAMPLE}_filtered.vcf.gz

        if [ -f "$FINAL_VCF" ]; then
            echo "$SAMPLE already complete → skipping"
            continue
        fi

        RECAL_BAM=$RESULTS/recal/${SAMPLE}_recal.bam

        if [ ! -f "$RECAL_BAM" ]; then
            echo "ERROR: $RECAL_BAM not found, skipping $SAMPLE"
            continue
        fi

        UNFILTERED=$RESULTS/mutect2/${SAMPLE}_unfiltered.vcf.gz
        TMP_VCF=$RESULTS/mutect2/${SAMPLE}_unfiltered.vcf.gz.tmp
        F1R2=$RESULTS/mutect2/${SAMPLE}_f1r2.tar.gz
        PRIORS=$RESULTS/mutect2/${SAMPLE}_artifact-priors.tar.gz
        PILEUP=$RESULTS/mutect2/${SAMPLE}_pileups.table
        CONTAM=$RESULTS/mutect2/${SAMPLE}_contamination.table

        echo "Running Mutect2 (tumor-only) for $SAMPLE"

        ################################
        # Mutect2 (tumor-only)
        ################################

        gatk --java-options "-Xmx4G" Mutect2 \
            -R $REF \
            -I $RECAL_BAM \
            -L $INTERVALS \
            -tumor $SAMPLE \
            --germline-resource $GNOMAD \
            --f1r2-tar-gz $F1R2 \
            -O $TMP_VCF \
            2> $LOGS/${SAMPLE}_mutect2.log

        if [ $? -eq 0 ]; then
            mv $TMP_VCF $UNFILTERED
        else
            echo "Mutect2 failed for $SAMPLE"
            rm -f $TMP_VCF
            continue
        fi

        ################################
        # Learn orientation artifacts
        ################################

        gatk --java-options "-Xmx4G" LearnReadOrientationModel \
            -I $F1R2 \
            -O $PRIORS

        ################################
        # Contamination estimation
        ################################
        
        gatk --java-options "-Xmx4G" GetPileupSummaries \
            -I $RECAL_BAM \
            -V $GNOMAD_COMMON \
            -L $INTERVALS \
            -O $PILEUP

        gatk --java-options "-Xmx4G" CalculateContamination \
            -I $PILEUP \
            -O $CONTAM

        ################################
        # Variant filtering
        ################################

        gatk --java-options "-Xmx4G" FilterMutectCalls \
            -R $REF \
            -V $UNFILTERED \
            --contamination-table $CONTAM \
            --ob-priors $PRIORS \
            -O $FINAL_VCF

        if [ $? -ne 0 ]; then
            echo "Filtering failed for $SAMPLE"
            continue
        fi

        ################################
        # CLEANUP TEMP FILES
        ################################
        rm -f $F1R2 $PRIORS $PILEUP $CONTAM
        echo "Temporary files for $SAMPLE removed"

        echo "$SAMPLE finished"

    done

    echo ""
    read -p "Batch finished. Press ENTER for next batch..."

done

echo "Pipeline finished"
