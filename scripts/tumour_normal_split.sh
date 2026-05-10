#!/bin/bash

INPUT="runinfo.csv"

TUMOR="tumor_samples.txt"
NORMAL="normal_samples.txt"
UNKNOWN="unknown_samples.txt"

# Clear outputs
> $TUMOR
> $NORMAL
> $UNKNOWN

awk -F',' '
NR==1 {
    # Find column indices dynamically
    for (i=1; i<=NF; i++) {
        if ($i=="Run") run_col=i
        if ($i=="Tumor") tumor_col=i
        if ($i=="Affection_Status") aff_col=i
        if ($i=="Disease") dis_col=i
    }
    next
}

{
    run=$run_col
    tumor=tolower($tumor_col)
    aff=tolower($aff_col)
    disease=tolower($dis_col)

    if (tumor ~ /yes|tumou?r|1/) {
        print run >> "'"$TUMOR"'"
    }
    else if (tumor ~ /no|normal|0/) {
        print run >> "'"$NORMAL"'"
    }
    else if (aff ~ /affected/) {
        print run >> "'"$TUMOR"'"
    }
    else if (aff ~ /unaffected/) {
        print run >> "'"$NORMAL"'"
    }
    else if (disease ~ /cancer|carcinoma|tumou?r/) {
        print run >> "'"$TUMOR"'"
    }
    else {
        print run >> "'"$UNKNOWN"'"
    }
}
' $INPUT

echo "Tumor: $(wc -l < $TUMOR)"
echo "Normal: $(wc -l < $NORMAL)"
echo "Unknown: $(wc -l < $UNKNOWN)"
