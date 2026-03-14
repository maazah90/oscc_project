#!/bin/bash

################################
# USER PATHS
################################
# Directory of this script
SCRIPT_DIR=$(dirname "$(realpath "$0")")
# Project root is assumed to be parent of scripts folder
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
mkdir -p "$RESULTS/pon"
mkdir -p "$RESULTS/mutect2"
mkdir -p "$LOGS"

################################
# SAMPLE LIST
################################
SAMPLES=($(cat $PROJECT_DIR/samples.txt))

TOTAL=${#SAMPLES[@]}
BATCH=4

################################
# LOOP OVER BATCHES
################################
for ((i=0;i<$TOTAL;i+=BATCH))
do
    echo "===================================="
    echo "Processing batch: ${SAMPLES[@]:i:BATCH}"
    echo "===================================="

    for SAMPLE in "${SAMPLES[@]:i:BATCH}"
    do
        RECAL_BAM=$RESULTS/recal/${SAMPLE}_recal.bam
        PON_VCF=$RESULTS/pon/${SAMPLE}_pon.vcf.gz
        TMP_VCF=$RESULTS/pon/${SAMPLE}_pon.vcf.gz.tmp

        # Skip finished samples
        if [ -f "$PON_VCF" ]; then
            echo "$SAMPLE already finished → skipping"
            continue
        fi

        # Check BAM exists
        if [ ! -f "$BAMDIR/${SAMPLE}.dedup.bam" ]; then
            echo "Error: $BAMDIR/${SAMPLE}.dedup.bam not found, skipping $SAMPLE"
            continue
        fi

        echo "Running $SAMPLE"

        ################################
        # BQSR
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


        gatk --java-options "-Xmx4G" ApplyBQSR \
            -R $REF \
            -I $BAMDIR/${SAMPLE}.dedup.bam \
            --bqsr-recal-file $RESULTS/recal/${SAMPLE}_recal.table \
            -O $RECAL_BAM \
            2>> $LOGS/${SAMPLE}_bqsr.log
            
          #Stop if BQSR failed  
            if [ $? -ne 0 ]; then
    		echo "BQSR failed for $SAMPLE"
    		continue
   	    fi


        ################################
        # Index recalibrated BAM
        ################################
        samtools index $RECAL_BAM

        ################################
        # Mutect2 PoN mode
        ################################
	gatk --java-options "-Xmx4G" Mutect2 \
	     -R $REF \
	     -I $RECAL_BAM \
	     -L $INTERVALS \
	     --max-mnp-distance 0 \
	     -O $TMP_VCF \
	     2> $LOGS/${SAMPLE}_pon.log



        # Only rename file if Mutect2 finished successfully
        if [ $? -eq 0 ]; then
            mv $TMP_VCF $PON_VCF
            echo "$SAMPLE finished successfully"
        else
            echo "$SAMPLE Mutect2 failed"
            rm -f $TMP_VCF
        fi

    done

    echo ""
    echo "Batch completed."
    read -p "Press ENTER to run the next batch..."

done

echo "Step 1 finished"

