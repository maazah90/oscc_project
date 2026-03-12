#!/bin/bash

################################
# USER PATHS
################################

REF=reference/GRCh38.fa
GNOMAD=resources/af-only-gnomad.vcf.gz
GNOMAD_COMMON=resources/gnomad_common.vcf.gz

INTERVALS=intervals/exome_targets.bed

RESULTS=results

mkdir -p $RESULTS/mutect2
mkdir -p logs

################################
# CREATE PANEL OF NORMALS
################################

if [ ! -f "$RESULTS/pon.vcf.gz" ]; then

echo "Creating Panel of Normals"

gatk --java-options "-Xmx8G" GenomicsDBImport \
--genomicsdb-workspace-path $RESULTS/pon_db \
$(for f in $RESULTS/pon/*.vcf.gz; do echo -V $f; done)

gatk --java-options "-Xmx8G" CreateSomaticPanelOfNormals \
-R $REF \
-V gendb://$RESULTS/pon_db \
-O $RESULTS/pon.vcf.gz

echo "PoN created"

fi

################################
# SAMPLE LIST
################################

SAMPLES=($(cat samples.txt))
TOTAL=${#SAMPLES[@]}
BATCH=4

################################
# FINAL MUTECT2 LOOP
################################

for ((i=0;i<$TOTAL;i+=BATCH))
do

echo "===================================="
echo "Running batch ${SAMPLES[@]:i:BATCH}"
echo "===================================="

for SAMPLE in "${SAMPLES[@]:i:BATCH}"
do

FINAL_VCF=$RESULTS/mutect2/${SAMPLE}_filtered.vcf.gz

if [ -f "$FINAL_VCF" ]; then
    echo "$SAMPLE already complete → skipping"
    continue
fi

RECAL_BAM=$RESULTS/recal/${SAMPLE}_recal.bam

################################
# Mutect2
################################

gatk --java-options "-Xmx4G" Mutect2 \
-R $REF \
-I $RECAL_BAM \
-L $INTERVALS \
-tumor $SAMPLE \
--panel-of-normals $RESULTS/pon.vcf.gz \
--germline-resource $GNOMAD \
-O $RESULTS/mutect2/${SAMPLE}_unfiltered.vcf.gz \
2> logs/${SAMPLE}_mutect2.log

################################
# contamination
################################

gatk --java-options "-Xmx4G" GetPileupSummaries \
-I $RECAL_BAM \
-V $GNOMAD_COMMON \
-L $INTERVALS \
-O $RESULTS/mutect2/${SAMPLE}_pileups.table

gatk --java-options "-Xmx4G" CalculateContamination \
-I $RESULTS/mutect2/${SAMPLE}_pileups.table \
-O $RESULTS/mutect2/${SAMPLE}_contamination.table

################################
# filtering
################################

gatk --java-options "-Xmx4G" FilterMutectCalls \
-R $REF \
-V $RESULTS/mutect2/${SAMPLE}_unfiltered.vcf.gz \
--contamination-table $RESULTS/mutect2/${SAMPLE}_contamination.table \
-O $FINAL_VCF

echo "$SAMPLE finished"

done

read -p "Batch finished. Press ENTER for next batch..."

done

echo "Pipeline finished"
