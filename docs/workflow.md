# OSCC Whole Exome Sequencing Analysis Pipeline
Author: [MAAZAH MUHAMMAD ALI]
Project: OSCC – Pakistani Cohort (Khyber Pakhtunkhwa)
Reference Build: GRCh38 (Ensembl primary assembly)
Platform: Ubuntu 22.04
Sequencing: Illumina WES
Samples: 40 tumour samples (paired-end FASTQ)

---

# 1. Project Structure

oscc_project/
│
├── fastq/               # Raw FASTQ files (.fastq.gz)
├── sra/                 # Downloaded SRA files
├── reference/           # GRCh38 reference genome
├── aligned/             # BAM files (sorted, deduplicated)
├── variants/            # VCF files
├── tmp/                 # Temporary files
├── logs/                # Log files
└── workflow.md          # This file

---

# 2. Data Download (SRA → FASTQ)

## Download SRA 
wget or prefetch SRRXXXXXXX.sra

## Convert to FASTQ
fasterq-dump SRRXXXXXXX.sra --split-files --threads 3 --outdir fastq/
gzip fastq/SRRXXXXXXX_*.fastq

## Download FASTQ directly

- Attempted to downloaded SRA and its dependencies from SRA NCBI but it kept failing. 
- Opted to download FastQ files from ENA which downloaded. The conversion to FastQ step was skipped over
---

# 3. Quality Control

## FastQC
fastqc fastq/*.fastq.gz -o qc/

## MultiQC
multiqc qc/

### Decision:
- Per-base quality: PASS
- Adapter content: minimal
- No trimming performed

---

# 4. Reference Genome Setup

## Download GRCh38 (Ensembl)

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.fa

## Indexing
bwa index GRCh38.fa
samtools faidx GRCh38.fa
picard CreateSequenceDictionary R=GRCh38.fa O=GRCh38.dict

---

# 5. Alignment (BWA MEM)

## Script: align_wes.sh

bwa mem -t 4 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"
GRCh38.fa
${SAMPLE}_1.fastq.gz
${SAMPLE}_2.fastq.gz
| samtools view -bS - > ${SAMPLE}.bam

samtools sort -@ 4 -o ${SAMPLE}.sorted.bam ${SAMPLE}.bam

picard MarkDuplicates
I=${SAMPLE}.sorted.bam
O=${SAMPLE}.dedup.bam
M=${SAMPLE}.metrics.txt

samtools index ${SAMPLE}.dedup.bam


## Parallel Execution (5 at a time)

cat samples.txt | parallel -j 5 bash align_wes.sh {}

System RAM: ~20GB  
Parallel jobs: 5 max
---

# 6. Base Quality Score Recalibration (BQSR)

Required resources:
- dbSNP (hg38)
- Mills and 1000G indels (hg38)

## BaseRecalibrator
gatk BaseRecalibrator
-R GRCh38.fa
-I ${SAMPLE}.dedup.bam
--known-sites dbsnp.vcf.gz
--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
-O ${SAMPLE}.recal.table

mkdir -p recal

for bam in aligned/markdup/*.markdup.bam
do
    base=$(basename $bam .markdup.bam)

    gatk BaseRecalibrator \
    -R reference/GRCh38.fa \
    -I $bam \
    --known-sites reference/dbsnp.vcf.gz \
    -O recal/${base}.recal.table

    gatk ApplyBQSR \
    -R reference/GRCh38.fa \
    -I $bam \
    --bqsr-recal-file recal/${base}.recal.table \
    -O recal/${base}.recal.bam

    samtools index recal/${base}.recal.bam
done

## ApplyBQSR
gatk ApplyBQSR
-R GRCh38.fa
-I ${SAMPLE}.dedup.bam
--bqsr-recal-file ${SAMPLE}.recal.table
-O ${SAMPLE}.recal.bam
---

# 7. Somatic Variant Calling (Tumour-Only Mode)

## Mutect2
gatk Mutect2
-R GRCh38.fa
-I ${SAMPLE}.recal.bam
-tumor ${SAMPLE}
--germline-resource gnomad.hg38.vcf.gz
-O ${SAMPLE}.unfiltered.vcf.gz

mkdir -p variants/raw

for bam in recal/*.recal.bam
do
    base=$(basename $bam .recal.bam)

    gatk Mutect2 \
    -R reference/GRCh38.fa \
    -I $bam \
    --germline-resource reference/gnomad.vcf.gz \
    -L reference/Agilent_V7_targets.bed \
    -O variants/raw/${base}.unfiltered.vcf.gz
done

## Filter Calls
gatk FilterMutectCalls
-V ${SAMPLE}.unfiltered.vcf.gz
-O ${SAMPLE}.filtered.vcf.gz
---

# 8. Variant Filtering Criteria

Minimum thresholds:
- Depth (DP) ≥ 20
- Genotype Quality (GQ) ≥ 20
- Quality Score (QUAL) ≥ 50

mkdir -p variants/filtered

for vcf in variants/raw/*.unfiltered.vcf.gz
do
    base=$(basename $vcf .unfiltered.vcf.gz)

    gatk FilterMutectCalls \
    -R reference/GRCh38.fa \
    -V $vcf \
    -O variants/filtered/${base}.filtered.vcf.gz
done

---

# 9. Annotation

Tool: ANNOVAR
table_annovar.pl sample.vcf humandb/
-buildver hg38
-out sample
-remove
-protocol refGene,clinvar_20220320,gnomad30_genome
-operation g,f,f
-nastring .
-csvout

mkdir -p annotated

for vcf in variants/filtered/*.filtered.vcf.gz
do
    base=$(basename $vcf .filtered.vcf.gz)

    vep \
    -i $vcf \
    -o annotated/${base}.vep.vcf \
    --cache \
    --assembly GRCh38 \
    --vcf \
    --symbol \
    --canonical \
    --protein \
    --af
done

Output: CSV for downstream analysis in R

## Merge all Samples

bcftools merge variants/filtered/*.vcf.gz -Oz -o variants/cohort_merged.vcf.gz
bcftools index variants/cohort_merged.vcf.gz
---

# 10. Gene Targets for Analysis

## Final Curated Gene Panel

1️⃣  WNT / β-Catenin Axis Canonical WNT Signalling

CTNNB1

APC

AXIN1

AXIN2

GSK3B

TCF7L2

### WNT Ligands & Receptors

WNT1

WNT3A

WNT5A

FZD1

FZD7

LRP6

### Downstream Targets / Proliferation

MYC

CCND1

Adhesion Crosstalk

CDH1

2️⃣ Core OSCC / HNSCC Driver Genes

TP53

NOTCH1

PIK3CA

FAT1

CDKN2A

HRAS

EGFR

CASP8

- Note: These are frequently altered in head & neck cancers.

3️⃣ DNA Damage / Genome Stability (For TMB Context)

- Important for tobacco-associated cancers.

ATM

ATR

BRCA1

BRCA2

CHEK2

4️⃣ PI3K/AKT/mTOR Pathway

Commonly altered in OSCC.

PIK3CA

PTEN

AKT1

MTOR

5️⃣ Additional High-Relevance OSCC Genes

FBXW7

NFE2L2

KEAP1

SMAD4

TGFBR2
---

# 11. Comparative Genomics

Potential comparisons:
- 1000 Genomes Punjabi (Lahore, Pakistan) from the 1000 Genomes Project
- European OSCC datasets from The Cancer Genome Atlas or COSMIC
- gnomAD population frequencies

Objectives:
- Identify population-enriched SNPs
- Compare somatic mutation profiles
- Evaluate mutational burden differences

## Compute allelic Frequency
bcftools isec -p comparative/pjl_overlap \
variants/cohort_merged.vcf.gz pjl.vcf.gz

---

# 12. Machine Learning (Exploratory)

Planned analyses:
- PCA on mutation matrix
- Logistic regression (mutation presence vs group)
- Mutation burden clustering
- Pathway enrichment modeling

Note: Small sample size → exploratory only

---

# 13. Reproducibility Notes

- Ubuntu 22.04
- BWA 
- Samtools
- Picard
- GATK
- GNU Parallel
- FastQC
- MultiQC
- ANNOVAR
- R / Python (ML phase)

All scripts version-controlled via Git.

## Exact Software Versions

BWA v0.7.17-r1188
samtools v1.13
Picard v3.4.0
Java version 21.0.10
Linux OS: Ubuntu 22.04 LTS
---

# 14. Publication Goal

Objective:
Population-aware somatic mutation profiling of OSCC in Pakistani cohort with WNT pathway focus.

Future expansion:
- Integrate normal controls
- Compare to global OSCC datasets
- Functional pathway enrichment
- Exploratory ML modeling

