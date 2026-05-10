#!/usr/bin/env python3
import os
from pathlib import Path

# === CONFIG ===
PLAIN_VCF_DIR = "/media/maazah/Expansion/oscc_project/results/plain_vcfs"
OUTDIR = "/media/maazah/Expansion/oscc_project/SigProfiler_results_GRCh38_WES"
MIN_SIG = 1
MAX_SIG = 3
NMF_REPLICATES = 100
EXOME = True
GENOME_NAME = "GRCh38"

SIGPROFILER_DATA_DIR = os.path.expanduser(
    "~/.local/lib/python3.10/site-packages/SigProfilerAssignment/data/Reference_Signatures"
)
GENOME_PATH = os.path.join(SIGPROFILER_DATA_DIR, GENOME_NAME)


def get_valid_vcfs(input_folder):
    """Return list of non-empty VCF files in input_folder."""
    input_path = Path(input_folder)
    if not input_path.is_dir():
        raise FileNotFoundError(f"Input VCF folder not found: {input_folder}")

    vcf_files = sorted(input_path.glob("*.vcf"))
    valid_vcfs = []
    for vcf in vcf_files:
        with open(vcf) as f:
            variant_count = sum(1 for line in f if not line.startswith("#") and line.strip())
        if variant_count > 0:
            valid_vcfs.append(str(vcf))
        else:
            print(f"Skipping empty VCF: {vcf.name}")
    return valid_vcfs


if __name__ == '__main__':
    from SigProfilerExtractor import sigpro as sp

    # Check genome exists
    if not os.path.isdir(GENOME_PATH):
        raise FileNotFoundError(f"{GENOME_NAME} genome not found at {GENOME_PATH}")

    # Get valid VCFs
    valid_vcfs = get_valid_vcfs(PLAIN_VCF_DIR)
    if not valid_vcfs:
        raise ValueError("No valid VCFs with variants found in plain VCF folder.")

    print(f"Using {len(valid_vcfs)} VCFs for SigProfilerExtractor.")
    print("Running SigProfilerExtractor on CPU...")

    sp.sigProfilerExtractor(
        "vcf",
        OUTDIR,
        PLAIN_VCF_DIR,
        reference_genome=GENOME_NAME,
        minimum_signatures=MIN_SIG,
        maximum_signatures=MAX_SIG,
        nmf_replicates=NMF_REPLICATES,
        exome=EXOME,
        cpu=True
    )

    print(f"Done. Results saved in: {OUTDIR}")
