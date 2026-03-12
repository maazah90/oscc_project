#!/bin/bash

################################

# USER PATHS

################################

REF=reference/GRCh38.fa
DBSNP=known_sites/dbsnp.vcf.gz
MILLS=known_sites/Mills_and_1000G_gold_standard.indels.vcf.gz

BAMDIR=aligned
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

echo "Running $SAMPLE"

################################

# BQSR

################################

gatk --java-options "-Xmx4G" BaseRecalibrator 
-R $REF 
-I $BAMDIR/${SAMPLE}.dedup.bam 
--known-sites $DBSNP 
--known-sites $MILLS 
-O $RESULTS/recal/${SAMPLE}_recal.table 
2> logs/${SAMPLE}_bqsr.log

gatk --java-options "-Xmx4G" ApplyBQSR 
-R $REF 
-I $BAMDIR/${SAMPLE}.dedup.bam 
--bqsr-recal-file $RESULTS/recal/${SAMPLE}_recal.table 
-O $RECAL_BAM 
2>> logs/${SAMPLE}_bqsr.log

################################

# Index recalibrated BAM

################################

samtools index $RECAL_BAM

################################

# Mutect2 PoN mode

################################

gatk --java-options "-Xmx4G" Mutect2 
-R $REF 
-I $RECAL_BAM 
--max-mnp-distance 0 
-O $TMP_VCF 
2> logs/${SAMPLE}_pon.log

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

