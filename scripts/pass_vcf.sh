#!/bin/bash

set -e
set -o pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

INPUT_DIR="$PROJECT_DIR/results/mutect2_trimmed"
OUTPUT_DIR="$PROJECT_DIR/results/pass_vcfs"

mkdir -p "$OUTPUT_DIR"

export OUTPUT_DIR

ls "$INPUT_DIR"/*_filtered.vcf.gz | xargs -n 1 -P 4 -I {} bash -c '
    FILE="{}"
    SAMPLE=$(basename "$FILE" _filtered.vcf.gz)

    echo "Processing $SAMPLE..."

    bcftools view -f PASS "$FILE" -Oz -o "$OUTPUT_DIR/${SAMPLE}_pass.vcf.gz"

    bcftools index "$OUTPUT_DIR/${SAMPLE}_pass.vcf.gz"

    echo "$SAMPLE done"
'

echo "All PASS VCFs generated"
