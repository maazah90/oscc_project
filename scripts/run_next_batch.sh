#!/bin/bash

THREADS=2
JOBS=2

BASE="/media/maazah/Expansion/oscc_project"
LOGDIR="$BASE/results/logs_trimmed"

mkdir -p "$LOGDIR"

if [ ! -f samples.txt ]; then
    echo "samples.txt not found"
    exit 1
fi

TOTAL=$(wc -l < samples.txt)
START=1

while [ $START -le $TOTAL ]; do

    echo "========================================="
    echo "Processing samples $START to $((START + JOBS - 1))"
    echo "========================================="

    # Extract current batch
    sed -n "${START},$((START + JOBS - 1))p" samples.txt > current_batch.txt

    if [ ! -s current_batch.txt ]; then
        echo "No more samples to process."
        break
    fi

    cat current_batch.txt
    echo "-----------------------------------------"

    parallel -j "$JOBS" --joblog "$LOGDIR/parallel_joblog.txt" \
    'echo "Starting {1} at $(date)";
     bash align_wes.sh {1} '"$THREADS"' > '"$LOGDIR"'/{1}.log 2>&1;
     echo "Finished {1} at $(date)"' \
    :::: current_batch.txt

    echo "-----------------------------------------"
    echo "Batch finished."

    # Ask for confirmation
    read -p "Run next batch? (y/n): " choice

    if [[ "$choice" != "y" ]]; then
        echo "Stopping pipeline."
        break
    fi

    START=$((START + JOBS))

done

echo "========================================="
echo "Pipeline ended."
