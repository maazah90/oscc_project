#!/bin/bash

BATCH_SIZE=2

head -n $BATCH_SIZE remaining.txt > current_batch.txt

if [ ! -s current_batch.txt ]; then
    echo "No more samples to process."
    exit 0
fi

echo "Running next batch:"
cat current_batch.txt
echo "-------------------------"

parallel -j 2 './align_wes.sh {1} 2' :::: current_batch.txt

# If successful, remove completed samples
tail -n +$(($BATCH_SIZE + 1)) remaining.txt > tmp.txt
mv tmp.txt remaining.txt

echo "Batch complete."
