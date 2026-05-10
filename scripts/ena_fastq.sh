#!/bin/bash
set -euo pipefail

ACCESSION="PRJEB14203"
OUTDIR="/media/maazah/Expansion/oscc_project/indian_study_samples"
THREADS="${3:-8}"

mkdir -p "$OUTDIR"

RAW="${OUTDIR}/ena_raw.tsv"

echo "🔎 Fetching ENA data..."

################################
# FETCH (DO NOT OVER-FILTER)
################################
curl -fsS \
"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
-o "$RAW"

################################
# VALIDATE DATA EXISTS
################################
TOTAL=$(tail -n +2 "$RAW" | wc -l)

if [[ "$TOTAL" -eq 0 ]]; then
    echo "❌ ERROR: No runs returned from ENA"
    exit 1
fi

echo "✔ Found $TOTAL runs"

################################
# DOWNLOAD FUNCTION
################################
download_run () {
    local RUN="$1"
    local LINKS="$2"

    local DIR="${OUTDIR}/${RUN}"
    mkdir -p "$DIR"

    # skip if already complete
    if compgen -G "${DIR}/*_R1.fastq.gz" >/dev/null && \
       compgen -G "${DIR}/*_R2.fastq.gz" >/dev/null; then
        echo "⏭️ Skipping $RUN"
        return
    fi

    echo "⬇️ Downloading $RUN"

    for url in $(echo "$LINKS" | tr ';' ' '); do
        wget -c -q --show-progress -P "$DIR" "$url"
    done

    cd "$DIR"

    for f in *_1.fastq.gz; do
        [[ -e "$f" ]] && mv -n "$f" "${RUN}_R1.fastq.gz"
    done

    for f in *_2.fastq.gz; do
        [[ -e "$f" ]] && mv -n "$f" "${RUN}_R2.fastq.gz"
    done

    cd - >/dev/null
}

################################
# PROCESS ALL RUNS (NO FILTERING)
################################
echo "🚀 Starting downloads..."

tail -n +2 "$RAW" | while IFS=$'\t' read -r RUN LINKS; do
    [[ -z "$RUN" || -z "$LINKS" ]] && continue
    download_run "$RUN" "$LINKS"
done

echo "✅ Done"
