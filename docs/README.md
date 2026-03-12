# OSCC Whole Exome Somatic Variant Analysis Pipeline

## Overview

This repository contains a computational workflow for detecting and analyzing somatic variants from Whole Exome Sequencing (WES) data derived from Oral Squamous Cell Carcinoma (OSCC) tumour samples.

The analysis includes alignment, somatic variant calling, filtering, functional annotation, and cohort-level mutation analysis.

A secondary objective is to explore potential population-specific mutation patterns by comparing the cohort with publicly available cancer datasets and population genomic resources.

---

# Study Objectives

* Identify somatic variants in OSCC tumour samples using a tumour-only pipeline.
* Determine genes frequently mutated across the cohort.
* Perform functional annotation of variants.
* Compare mutation patterns with publicly available cancer cohorts.
* Explore potential population-specific genomic signatures.

---

# Dataset

Tumour-only Whole Exome Sequencing samples.

Number of samples: **40**

Sample IDs:

```text
SRR31443517 – SRR31443556
```

Reference genome: **GRCh38**

Data source: Sequence Read Archive (SRA)

---

# Analysis Pipeline

## Pipeline Overview

```text
Raw FASTQ
   │
   ▼
Alignment (BWA)
   │
   ▼
Sort + MarkDuplicates
   │
   ▼
BQSR
   │
   ▼
Mutect2 Variant Calling
   │
   ▼
FilterMutectCalls
   │
   ▼
Optional Variant Filtering
   │
   ▼
ANNOVAR Annotation
   │
   ▼
Gene-Level Mutation Analysis
   │
   ▼
Cohort Comparison
   │
   ▼
Optional Population PCA
```

---

# Workflow Description

## 1. Read Alignment

Sequencing reads are aligned to the human reference genome using **BWA**.

Post-alignment processing includes:

* sorting reads
* marking PCR duplicates
* indexing BAM files

Tools used:

* **SAMtools**
* **Picard**

Output:

```text
sample_dedup.bam
sample_dedup.bam.bai
```

---

# 2. Base Quality Score Recalibration

Sequencing error patterns are corrected using **GATK** Base Quality Score Recalibration.

Steps:

1. BaseRecalibrator
2. ApplyBQSR

Known variant databases used:

* dbSNP
* Mills and 1000G indels

---

# 3. Somatic Variant Calling

Somatic mutation detection is performed using **Mutect2**.

Pipeline steps:

1. Mutect2
2. FilterMutectCalls

Because matched normal samples are unavailable, a tumour-only analysis approach is used.

A **Panel of Normals (PoN)** is generated to remove recurrent technical artefacts.

Output:

```text
results/mutect2/
sample_filtered.vcf.gz
```

---

# 4. Variant Filtering

Additional variant filtering may be applied depending on study requirements.

Example thresholds:

| Metric                | Threshold |
| --------------------- | --------- |
| Depth (DP)            | ≥ 20      |
| Quality Score (QUAL)  | ≥ 50      |
| Genotype Quality (GQ) | ≥ 20      |

Filtering can be performed using **bcftools**.

---

# 5. Variant Annotation

Functional annotation is performed using **ANNOVAR**.

Annotation includes:

* affected genes
* variant functional effect
* population allele frequency
* clinical significance

Databases commonly used:

* RefSeq
* ClinVar
* gnomAD
* dbSNP

---

# 6. Gene-Level Mutation Analysis

Annotated variants are analyzed to identify genes frequently mutated across samples.

Common OSCC-associated genes include:

* TP53
* NOTCH1
* FAT1
* PIK3CA
* CDKN2A

Possible analyses include:

* mutation frequency per gene
* mutation type distribution
* recurrent mutation identification

---

# 7. Cohort-Level Mutation Analysis

After variant annotation and filtering, cohort-level mutation analysis will be performed using **maftools** in **R**.

Annotated variants will be converted into Mutation Annotation Format (MAF) files, which are commonly used for cancer genomics studies.

This analysis enables visualization and statistical exploration of somatic mutation patterns across the cohort.

### Analyses performed with maftools

Possible analyses include:

* Mutation frequency analysis
* Gene mutation summaries
* Mutation spectrum analysis
* Tumour mutation burden estimation
* Cohort visualization

### Common visualizations

Typical visualizations generated using maftools include:

* Oncoplots (mutation landscape across samples)
* Mutation frequency barplots
* Variant classification summaries
* Gene mutation heatmaps

These analyses help identify recurrently mutated genes and characterize the overall mutation profile of the cohort.

### Example OSCC driver genes examined

Frequently mutated genes in **Oral Squamous Cell Carcinoma** and related to my work on the **Wnt-Beta Catenin pathway**:

* TP53
* NOTCH1
* FAT1
* PIK3CA
* CDKN2A
* B-Catenin
* E-Cadherin
* c-MYC
* Cyclin D1

Mutation patterns observed in the cohort will be compared with public datasets such as **TCGA Head and Neck Squamous Cell Carcinoma (TCGA-HNSC)** and **COSMIC**.



# 8. Cohort Comparison

Mutation profiles may be compared with publicly available cancer datasets including:

* **TCGA Head and Neck Squamous Cell Carcinoma (TCGA-HNSC)**
* **COSMIC**

Comparisons may examine:

* mutation frequencies
* shared driver mutations
* mutation spectrum differences

---

# 9. Population-Level Analysis (Optional)

To explore potential population-specific genomic variation, the cohort may be compared with global population datasets including the **1000 Genomes Project**.

Variants may be converted into genotype format and analyzed using **PLINK**.

Principal Component Analysis (PCA) can then be performed to examine clustering patterns relative to known population groups.

This exploratory analysis may provide insight into potential genomic variation in underrepresented South Asian populations.

---

# Software Requirements

Required software:

* **BWA**
* **SAMtools**
* **Picard**
* **GATK**
* **ANNOVAR**
* **bcftools**
* **PLINK**
* **R**

---

# Repository Structure

```text
reference/        reference genome
known_sites/      dbSNP and indel databases
scripts/          pipeline scripts
results/          output files
logs/             runtime logs
```

---

# Running the Pipeline

The variant calling workflow is executed using bash scripts stored in the `scripts/` directory.

Example execution:

```bash
bash run_mutect2_batches.sh
```

After variant calling is complete, the Panel of Normals is generated:

```bash
bash create_panel_of_normals.sh
```

Intermediate files are retained to allow recovery if the pipeline is interrupted.

---

# Hardware Environment

Analysis performed on a personal workstation with:

* 4 CPU cores
* ~20 GB RAM
* external storage for sequencing data

Batch processing was used to prevent resource exhaustion.

---

# Logging

Runtime logs were recorded during alignment and variant calling.

These files track:

* sample processing status
* runtime per sample
* pipeline errors

Logs are excluded from version control using `.gitignore`.

---

# Future Directions

Potential future analyses include:

* mutation signature analysis
* tumour mutation burden estimation
* oncoplot visualization
* pathway enrichment analysis
* expanded cohort comparison

---

# Author
MAAZAH MUHAMMAD ALI
Computational genomics analysis of somatic mutations in OSCC Whole Exome Sequencing data.

