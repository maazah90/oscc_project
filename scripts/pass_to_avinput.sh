#!/bin/bash

set -euo pipefail

################################
# PATHS
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

VCF_DIR="$PROJECT_DIR/results/pass_vcfs_indian"
ANNOVAR_BIN="$PROJECT_DIR/annovar"
OUT_DIR="$PROJECT_DIR/annovar/avinput_indian"
LOGS="$OUT_DIR/logs"

mkdir -p "$OUT_DIR" "$LOGS"

################################
# LOOP
################################

for vcf in "$VCF_DIR"/*.pass.vcf.gz; do

    [[ -e "$vcf" ]] || { echo "❌ No VCF files found"; break; }

    sample=$(basename "$vcf" .pass.vcf.gz)
    avinput="$OUT_DIR/${sample}.avinput"

    echo "======================================"
    echo "🔍 Processing: $sample"

    ################################
    # DIRECT CONVERSION (TUMOR-ONLY SAFE)
    ################################
    echo "⏳ Converting to avinput..."

    zcat "$vcf" | \
    perl "$ANNOVAR_BIN/convert2annovar.pl" \
        -format vcf4 \
        - \
        > "$avinput" 2> "$LOGS/${sample}.log"

    ################################
    # VALIDATION
    ################################
    if [[ -s "$avinput" ]]; then
        echo "✅ Done: $sample ($(wc -l < "$avinput") variants)"
    else
        echo "❌ Failed: $sample"
        echo "   Check log: $LOGS/${sample}.log"
    fi

done

echo "🎉 All samples processed!"
