#!/bin/bash

THREADS=2
JOBS=2

BASE="/media/maazah/Expansion/oscc_project"
LOGDIR="$BASE/results/logs_indian_samples"

mkdir -p "$LOGDIR"

if [ ! -f indian_srrs_clean.txt ]; then
    echo "indian_srrs_clean.txt not found"
    exit 1
fi

TOTAL=$(wc -l < indian_srrs_clean.txt)
START=1

while [ $START -le $TOTAL ]; do

    echo "========================================="
    echo "Processing samples $START to $((START + JOBS - 1))"
    echo "========================================="

    # Extract current batch
    sed -n "${START},$((START + JOBS - 1))p" indian_srrs_clean.txt > current_batch_indian.txt

    if [ ! -s current_batch_indian.txt ]; then
        echo "No more samples to process."
        break
    fi

    cat current_batch_indian.txt
    echo "-----------------------------------------"

    parallel -j "$JOBS" --joblog "$LOGDIR/parallel_joblog.txt" \
    'echo "Starting {1} at $(date)";
     bash align_wes.sh {1} '"$THREADS"' > '"$LOGDIR"'/{1}.log 2>&1;
     echo "Finished {1} at $(date)"' \
    :::: current_batch_indian.txt

    echo "-----------------------------------------"
    echo "Batch finished."


    START=$((START + JOBS))

done

echo "========================================="
echo "Pipeline ended."
