#!/bin/bash

################################
# TUMOR-ONLY PIPELINE WRAPPER
################################
# Usage: ./run_pipeline.sh
# Make sure script1.sh and script2.sh are in the same folder
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

SCRIPT1=$SCRIPT_DIR/indian_samples_bqsr.sh   # BQSR
SCRIPT2=$SCRIPT_DIR/indian_mutect2_run.sh   # Tumor-only Mutect2

echo "===================================="
echo "Tumor-only WES Pipeline Wrapper"
echo "Project dir: $PROJECT_DIR"
echo "===================================="

################################
# STEP 1: BQSR
################################
echo ""
echo "Step 1: Running BQSR (script1)"
echo "------------------------------------"
bash $SCRIPT1
echo "Step 1 completed (or paused batch-by-batch)."
echo "Check results/recal_indian/ for *_recal.bam files."

# read -p "Press ENTER to continue to Mutect2 tumor-only..."

################################
# STEP 2: Mutect2 tumor-only
################################
echo ""
echo "Step 2: Running Mutect2 tumor-only (script2)"
echo "------------------------------------"
bash $SCRIPT2
echo "Step 2 completed (or paused batch-by-batch)."
echo "Check results/mutect2_indian_samples/ for *_filtered.vcf.gz files."

echo ""
echo "Pipeline finished successfully!"
