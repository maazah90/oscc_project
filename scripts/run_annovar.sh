#!/bin/bash

set -euo pipefail

################################
# PATHS
################################

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$(dirname "$SCRIPT_DIR")

AVINPUT_DIR="$PROJECT_DIR/annovar/avinput_tumour"
HUMANDB="$PROJECT_DIR/annovar/humandb"
OUTPUT_DIR="$PROJECT_DIR/annovar_detailed"
LOGS="$PROJECT_DIR/annovar/logs_detailed"

mkdir -p "$OUTPUT_DIR" "$LOGS"

################################
# RUN ANNOVAR
################################

for av in "$AVINPUT_DIR"/*.avinput; do

    [[ -e "$av" ]] || { echo "❌ No avinput files found"; break; }

    sample=$(basename "$av" .avinput)

    if [[ -s "$OUTPUT_DIR/${sample}.hg38_multianno.txt" ]]; then
        echo "✅ Already annotated: $sample"
        continue
    fi

    echo "🧬 Annotating: $sample"

 perl "$PROJECT_DIR/annovar/table_annovar.pl" \
    "$av" \
    "$HUMANDB" \
    -buildver hg38 \
    -out "$OUTPUT_DIR/$sample" \
    -remove \
    -protocol ensGene,clinvar_20220320,gnomad211_exome,dbnsfp42a,1000g2015aug_all,1000g2015aug_sas,cytoBand,cosmic70 \
    -operation g,f,f,f,f,f,r,f \
    -nastring . \
    -polish \
    > "$LOGS/${sample}.log" 2>&1

    if [[ -s "$OUTPUT_DIR/${sample}.hg38_multianno.txt" ]]; then
        echo "✅ Done: $sample"
    else
        echo "❌ Failed: $sample (check log)"
    fi

done

echo "🎉 All ANNOVAR jobs completed"
