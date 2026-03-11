# OSCC WES Analysis Pipeline

## Overview
Whole Exome Sequencing analysis pipeline for OSCC samples.

## Reference Genome
GRCh38

## Pipeline Steps
1. BWA MEM alignment
2. SAM → BAM conversion
3. Sorting
4. MarkDuplicates (Picard)
5. BAM indexing
6. Variant calling (GATK - planned)

## Script Usage

Alignment:
bash scripts/align_wes.sh SAMPLE_ID

Example:
bash scripts/align_wes.sh SRR31443518

## Hardware
Laptop, ~20GB RAM

## Notes
Running samples sequentially or in small batches to avoid RAM overload.
