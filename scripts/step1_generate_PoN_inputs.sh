#!/bin/bash

################################
# USER PATHS
################################

REF=reference/GRCh38.fa
DBSNP=known_sites/dbsnp.vcf.gz
MILLS=known_sites/Mills_and_1000G_gold_standard.indels.vcf.gz

INTERVALS=intervals/exome_targets.bed

BAMDIR=bam
RESULTS=results

mkdir -p $RESULTS/recal
mkdir -p $RESULTS/pon
mkdir -p logs

################################
# SAMPLE LIST
################################

SAMPLES=($(cat samples.txt))
TOTAL=${#SAMPLES[@]}
BATCH=4

################################
# LOOP OVER BATCHES
################################

for ((i=0;i<$TOTAL;i+=BATCH))
do

echo "===================================="
echo "Processing batch ${SAMPLES[@]:i:BATCH}"
echo "===================================="

for SAMPLE in "${SAMPLES[@]:i:BATCH}"
do

RECAL_BAM=$RESULTS/recal/${SAMPLE}_recal.bam
PON_VCF=$RESULTS/pon/${SAMPLE}_pon.vcf.gz

# Skip finished samples
if [ -f "$PON_VCF" ]; then
    echo "$SAMPLE already finished → skipping"
    continue
fi

echo "Running $SAMPLE"

################################
# BQSR
################################

gatk --java-options "-Xmx4G" BaseRecalibrator \
-R $REF \
-I $BAMDIR/${SAMPLE}_dedup.bam \
--known-sites $DBSNP \
--known-sites $MILLS \
-L $INTERVALS \
-O $RESULTS/recal/${SAMPLE}_recal.table \
2> logs/${SAMPLE}_bqsr.log

gatk --java-options "-Xmx4G" ApplyBQSR \
-R $REF \
-I $BAMDIR/${SAMPLE}_dedup.bam \
--bqsr-recal-file $RESULTS/recal/${SAMPLE}_recal.table \
-L $INTERVALS \
-O $RECAL_BAM \
2>> logs/${SAMPLE}_bqsr.log

################################
# Mutect2 PoN mode
################################

gatk --java-options "-Xmx4G" Mutect2 \
-R $REF \
-I $RECAL_BAM \
-L $INTERVALS \
--max-mnp-distance 0 \
-O $PON_VCF \
2> logs/${SAMPLE}_pon.log

echo "$SAMPLE finished"

done

echo ""
read -p "Batch complete. Press ENTER to continue..."

done

echo "Step 1 finished"
