#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Use the already existing plain VCF folder
PLAIN_VCF_DIR="$PROJECT_DIR/results/plain_vcfs"
OUTDIR="$PROJECT_DIR/SigProfiler_results_GRCh38_WES"

MIN_SIG=1
MAX_SIG=8
NMF_REPLICATES=50
EXOME=true
GENOME_NAME="GRCh38"

export SIGPROFILER_DATA_DIR="${HOME}/.local/lib/python3.10/site-packages/SigProfilerAssignment/data/Reference_Signatures"

# Ensure the genome is installed
if [ ! -d "$SIGPROFILER_DATA_DIR/$GENOME_NAME" ]; then
    echo "Error: $GENOME_NAME genome not found in $SIGPROFILER_DATA_DIR" >&2
    exit 1
fi

# Check plain VCF folder exists
if [ ! -d "$PLAIN_VCF_DIR" ]; then
    echo "Error: Plain VCF directory not found: $PLAIN_VCF_DIR" >&2
    exit 1
fi

# Step 1: Run SigProfiler on plain VCFs
python3 - <<PYTHON
import os
from pathlib import Path
from SigProfilerExtractor import sigpro as sp

INPUT_DIR = os.path.abspath("$PLAIN_VCF_DIR")
OUTDIR = os.path.abspath("$OUTDIR")
MIN_SIG = int("$MIN_SIG")
MAX_SIG = int("$MAX_SIG")
NMF_REPLICATES = int("$NMF_REPLICATES")
EXOME = True if "$EXOME".lower() in ("1", "true", "yes") else False
GENOME_NAME = "$GENOME_NAME"

# Verify genome path
GENOME_PATH = os.path.join(os.environ["SIGPROFILER_DATA_DIR"], GENOME_NAME)
if not os.path.isdir(GENOME_PATH):
    raise FileNotFoundError(f"{GENOME_NAME} genome not found at {GENOME_PATH}")

print(f"Using {GENOME_NAME} genome at: {GENOME_PATH}")

# Verify input files
input_path = Path(INPUT_DIR)
vcf_files = sorted(input_path.glob("*.vcf"))
valid_vcfs = []

for vcf in vcf_files:
    with open(vcf) as f:
        variant_count = sum(1 for line in f if not line.startswith("#") and line.strip())
    if variant_count > 0:
        valid_vcfs.append(str(vcf))
    else:
        print(f"Skipping empty VCF: {vcf.name}")

if not valid_vcfs:
    raise ValueError("No valid VCFs with variants found in plain VCF folder.")

print(f"Using {len(valid_vcfs)} VCFs for SigProfilerExtractor.")

print("Running SigProfilerExtractor on CPU...")
sp.sigProfilerExtractor(
    "vcf",
    OUTDIR,
    INPUT_DIR,
    reference_genome=GENOME_NAME,
    minimum_signatures=MIN_SIG,
    maximum_signatures=MAX_SIG,
    nmf_replicates=NMF_REPLICATES,
    exome=EXOME,
    cpu=True
)
print(f"Done. Results saved in: {OUTDIR}")
PYTHON

echo "Script completed successfully."
