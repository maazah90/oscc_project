# oscc_project
Code for the Comparative Genomics Paper

# Pakistani OSCC Comparative Genomic Analysis

## Overview

This repository contains the bioinformatics pipeline scripts and R 
analysis code for the study:

**"Distinct Mutational Landscapes and Increased Actionable Alterations 
in Pakistani Oral Squamous Cell Carcinoma"**

Maazah Muhammad Ali, Ziauddin University, Karachi, Pakistan

> Preprint available at: [bioRxiv link — add when posted]

---

## Study Summary

This study performs a comparative whole-exome sequencing (WES) based 
genomic analysis of a Pakistani OSCC cohort against TCGA-derived 
European oral cavity samples, with independent validation in an Indian 
tongue squamous cell carcinoma cohort. Key findings include elevated 
tumour mutational burden, absence of APOBEC-associated signatures, 
enrichment of DNA damage response pathway alterations, and a 
significantly higher proportion of clinically actionable alterations 
in the Pakistani cohort.

---

## Repository Structure
---

## Data Sources

| Cohort | Accession | Description |
|---|---|---|
| Pakistani OSCC | PRJNA1189482 | 29 WES tumour samples |
| Indian validation | ERP015832 | 24 WES FFPE tongue SCC samples |
| TCGA European | GDC portal | TCGA-HNSC oral cavity subset |

---

## Requirements

### System
- Ubuntu 22.04 LTS
- 4+ CPU cores recommended
- 16GB+ RAM recommended

### Bioinformatics Tools
| Tool | Version |
|---|---|
| BWA | v0.7.17 |
| Samtools | v1.13 |
| Picard | v3.4.0 |
| GATK | Java v21.0.10 |
| FastQC | latest |
| ANNOVAR | latest |
| GNU Parallel | latest |

### R Packages
| Package | Version |
|---|---|
| R | v4.5.3+ |
| maftools | Bioconductor |
| TCGAbiolinks | Bioconductor |
| ggplot2 | CRAN |
| dplyr | CRAN |
| data.table | CRAN |
| uwot | CRAN |
| ggpubr | CRAN |

Install R packages with:
```r
install.packages("BiocManager")
BiocManager::install(c("maftools", "TCGAbiolinks"))
install.packages(c("ggplot2", "dplyr", "tidyr", 
                   "data.table", "uwot", "ggpubr",
                   "ggrepel", "patchwork", "stringr"))
```

---

## Usage

### 1. WES Processing Pipeline

See `workflow.md` for the complete step-by-step pipeline.

Run scripts in order for each cohort:

```bash
# Pakistani cohort
bash pipeline/pakistani_cohort/01_alignment.sh
bash pipeline/pakistani_cohort/02_markdup_bqsr.sh
bash pipeline/pakistani_cohort/03_mutect2.sh
bash pipeline/pakistani_cohort/04_annovar.sh

# Indian cohort (includes FFPE orientation bias correction)
bash pipeline/indian_cohort/01_alignment.sh
bash pipeline/indian_cohort/02_markdup_bqsr.sh
bash pipeline/indian_cohort/03_mutect2_ffpe.sh
bash pipeline/indian_cohort/04_annovar.sh
```

### 2. R Analysis

Set your working directory to the repository root and run scripts 
in order:

```r
source("R_analysis/01_main_pipeline.R")
source("R_analysis/02_drug_gene_targets.R")
source("R_analysis/03_pjl_comparison.R")
```

All figures are saved to `FIGURES/` directory.

---

## Reference Genome

- Pakistani and Indian cohorts: GRCh38/hg38
- TCGA data: GRCh37 (MC3 dataset)

Note: All comparative analyses were performed at the gene level 
and are not affected by this reference genome difference. 
See manuscript Limitations section for full discussion.

---

## Panel of Normals

Tumour-only variant calling used the GATK best practices 
panel of normals:

Available from: 
https://storage.googleapis.com/gatk-best-practices/somatic-hg38/

---

## Citation

If you use these scripts or data in your research, please cite:

> Ali, M.M. Distinct Mutational Landscapes and Increased Actionable 
> Alterations in Pakistani Oral Squamous Cell Carcinoma. 
> bioRxiv [year]. doi: [add when posted]

---

## Contact

**Maazah Muhammad Ali**  
College of Molecular Medicine  
Ziauddin University, Karachi, Pakistan  
GitHub: [@maazah90](https://github.com/maazah90)

---

## License

This project is available under the MIT License. 
See LICENSE file for details. 
