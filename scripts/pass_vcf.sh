#!/bin/bash

set -e
set -o pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

INPUT_DIR="$PROJECT_DIR/results/mutect2_paired"
OUTPUT_DIR="$PROJECT_DIR/results/pass_vcfs_paired"

mkdir -p "$OUTPUT_DIR"

export OUTPUT_DIR

ls "$INPUT_DIR"/*.filtered.vcf.gz | xargs -n 1 -P 4 -I {} bash -c '
    FILE="{}"
    SAMPLE=$(basename "$FILE" .filtered.vcf.gz)
    OUT_VCF="$OUTPUT_DIR/${SAMPLE}.pass.vcf.gz"

    # ✅ Skip if VCF + index (csi or tbi) exist
    if [[ -f "$OUT_VCF" && ( -f "${OUT_VCF}.csi" || -f "${OUT_VCF}.tbi" ) ]]; then
        echo "✅ Skipping $SAMPLE (already exists)"
        exit 0
    fi

    echo "⏳ Processing $SAMPLE..."

    bcftools view -f PASS "$FILE" -Oz -o "$OUT_VCF"

    bcftools index "$OUT_VCF"

    echo "✅ $SAMPLE done"
'

    echo "✅ $SAMPLE done"
'
