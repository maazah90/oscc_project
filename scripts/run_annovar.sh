#!/bin/bash

set -e
set -o pipefail

################################
# PATHS
################################

PROJECT_DIR="/media/maazah/Expansion/oscc_project"
AVINPUT_DIR="$PROJECT_DIR/annovar/avinput"
HUMANDB="$PROJECT_DIR/annovar/humandb"
OUTPUT_DIR="$PROJECT_DIR/annovar_results"

mkdir -p "$OUTPUT_DIR"

################################
# RUN ANNOVAR
################################

for av in "$AVINPUT_DIR"/*.avinput; do

    sample=$(basename "$av" .avinput)
    if [[ -f "$OUTPUT_DIR/${sample}.hg38_multianno.txt" ]]; then
        echo "✅ Already annotated: $sample"
        continue
    fi

    echo "🧬 Annotating: $sample"

perl "$AVINPUT_DIR/table_annovar.pl" "$av" "$HUMANDB" \
    -buildver hg38 \
    -out "$OUTPUT_DIR/$sample" \
    -remove \
    -protocol ensGene,clinvar_20220320,gnomad211_exome \
    -operation g,f,f \
    -nastring . \
    -polish \
    > "$OUTPUT_DIR/${sample}.log" 2>&1
      

done

echo "$SAMPLE done"

'

echo "All ANNOVAR jobs completed"
