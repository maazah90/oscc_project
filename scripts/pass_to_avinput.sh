#!/bin/bash

# ============================================
# 🗂 PATHS
# ============================================
VCF_DIR="/media/maazah/Expansion/oscc_project/results/pass_vcfs"       # folder with your pass.vcf files
ANNOVAR_DIR="/media/maazah/Expansion/oscc_project/annovar/avinput"      # folder to save avinput files
HUMANDB="/media/maazah/Expansion/oscc_project/annovar/humandb"

mkdir -p "$ANNOVAR_DIR"

# ============================================
# 🔁 LOOP THROUGH VCF FILES
# ============================================

for vcf in "$VCF_DIR"/*_pass.vcf.gz; do

    [[ -e "$vcf" ]] || { echo "❌ No VCF files found"; break; }

    sample=$(basename "$vcf" _pass.vcf.gz)
    avinput="$ANNOVAR_DIR/${sample}.avinput"

    if [[ -f "$avinput" ]]; then
        echo "✅ Already done: $sample"
        continue
    fi

    echo "⏳ Converting: $sample"

    zcat "$vcf" | perl "$ANNOVAR_DIR/convert2annovar.pl" \
        -format vcf4 \
        -includeinfo \
        - \
        > "$avinput" 2> "${avinput}.log"

    # Check if output is empty
    if [[ ! -s "$avinput" ]]; then
        echo "❌ Failed: $sample (check ${avinput}.log)"
    else
        echo "✅ Done: $sample"
    fi

done
