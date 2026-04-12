#!/usr/bin/env bash
set -euo pipefail

# --- Paths ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PLAIN_VCF_DIR="$PROJECT_DIR/results/plain_vcfs"
OUTDIR="$PROJECT_DIR/SigProfiler_results_GRCh38_WES"

# Ensure output directory exists
mkdir -p "$OUTDIR"

echo "=== Running SigProfiler ==="
echo "Using plain VCFs from: $PLAIN_VCF_DIR"
echo "Results will be saved in: $OUTDIR"

# Run the Python script
python3 "$SCRIPT_DIR/run_sigprofiler.py"

echo "=== SigProfiler workflow completed ==="
