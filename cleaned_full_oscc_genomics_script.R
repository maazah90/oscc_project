############################################################
# OSCC COMPARATIVE GENOMICS — PUBLICATION PIPELINE
# South Asian vs TCGA European OSCC (hg38)
############################################################

cat("\n===== STARTING PUBLICATION PIPELINE =====\n\n")

############################################################
# 0. LIBRARIES (safe install order)
############################################################
.libs_paths <- c("~/Rlibs", "/media/maazah/Expansion/R_libs", .Library)
.libPaths(.libs_paths)

safe_install <- function(pkg, bioc=FALSE, lib=NULL){
  if(!require(pkg, character.only=TRUE)){
    if(bioc) BiocManager::install(pkg, ask=FALSE, update=FALSE, lib=lib)
    else install.packages(pkg, lib=lib)
    library(pkg, character.only=TRUE)
  }
}

# Core packages
safe_install("maftools", bioc=TRUE, lib="~/Rlibs")
safe_install("TCGAbiolinks", bioc=TRUE, lib="~/Rlibs")
safe_install("dplyr")
safe_install("stringr")
safe_install("tidyr")
safe_install("ggplot2")
safe_install("data.table")
safe_install("BSgenome.Hsapiens.UCSC.hg38", bioc=TRUE)

# ============================================
# 📁 PATHS
# ============================================
annovar_dir <- "annovar_results"
output_dir <- "FIGURES"
dir.create(output_dir, showWarnings = FALSE)

files <- list.files(annovar_dir, pattern = "*.multianno.txt", full.names = TRUE)

# ============================================
# 🧬 PATHWAY GENE SETS
# ============================================
wnt_genes <- toupper(c("CTNNB1","APC","AXIN1","AXIN2","TCF7L2","LRP5","LRP6","DKK1"))
pi3k_genes <- toupper(c("PIK3CA","PIK3R1","AKT1","AKT2","PTEN","MTOR","TSC1","TSC2"))
ras_genes <- toupper(c("KRAS","NRAS","HRAS","BRAF","MAP2K1","MAP2K2","MAPK1"))


############################################################
# 1. DEFINE DRIVER AND NOISE GENES
############################################################
driver_genes <- toupper(c(
  "TP53","CDKN2A","NOTCH1","FAT1",
  "CTNNB1","APC","AXIN1","AXIN2","TCF7L2",
  "PIK3CA","PIK3R1","AKT1","PTEN","MTOR",
  "KRAS","NRAS","HRAS","BRAF"
))

noise_genes <- toupper(c(
  "TTN","MUC16","MUC4","MUC5AC","MUC5B",
  "OBSCN","DST","FLG","DNAH1","DNAH2",
  "DNAH5","DNAH8","RYR1","AHNAK2",
  "PCLO","UBR4","LRP1","HMCN2",
  "EPPK1","MUC12","PLEC","SPTBN5",
  "APOB","SBF2","BLTP1","CSMD2"
))

clean_gene <- function(x){
  x <- as.character(x)
  x <- sapply(strsplit(x, "[,;]"), `[`, 1)
  x <- str_trim(x)
  x <- toupper(x)
  x[x %in% c("", ".", "NA")] <- NA
  return(x)
}

############################################################
# 2. MAF FUNCTION WITH FILTERING
############################################################
convert_annovar_to_df <- function(file){
  
  df <- fread(file, data.table = FALSE)
  
  # Ensure required columns exist
  required_cols <- c("Chr","Start","End","Ref","Alt","Gene.ensGene","ExonicFunc.ensGene")
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0){
    stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))
  }
  
  # 🔥 FORCE CONSISTENT TYPES (CRITICAL FIX)
  df$Chr <- as.character(df$Chr)
  df$Start <- as.numeric(df$Start)
  df$End <- as.numeric(df$End)
  df$Ref <- as.character(df$Ref)
  df$Alt <- as.character(df$Alt)
  df$Gene.ensGene <- as.character(df$Gene.ensGene)
  df$ExonicFunc.ensGene <- as.character(df$ExonicFunc.ensGene)
  
  # Safe filtering (handles missing AF columns)
  df <- df %>%
    filter(
      ExonicFunc.ensGene != "synonymous SNV",
      (!"AF" %in% colnames(df) | is.na(AF) | AF < 0.01),
      (!"AF_nfe" %in% colnames(df) | is.na(AF_nfe) | AF_nfe < 0.01)
    )
  
  # Safe ClinVar handling
  df <- df %>%
    mutate(
      ClinVar_Pathogenic = ifelse(
        "CLNSIG" %in% colnames(df),
        grepl("Pathogenic", CLNSIG, ignore.case = TRUE),
        FALSE
      ),
      ClinVar_Significance = ifelse(
        "CLNSIG" %in% colnames(df),
        CLNSIG,
        NA
      )
    )
  
  # Sample name
  sample <- tools::file_path_sans_ext(basename(file))
  sample <- sub("\\.hg38_multianno$", "", sample)
  
  if(grepl("_vs_", sample)){
    sample <- strsplit(sample, "_vs_")[[1]][1]
  }
  
  # Build MAF-like dataframe
  maf <- data.frame(
    Hugo_Symbol = clean_gene(df$Gene.ensGene),
    Chromosome = as.character(df$Chr),   # 🔥 ALWAYS character
    Start_Position = as.numeric(df$Start),
    End_Position = as.numeric(df$End),
    Reference_Allele = toupper(df$Ref),
    Tumor_Seq_Allele2 = toupper(df$Alt),
    Variant_Classification = df$ExonicFunc.ensGene,
    Variant_Type = ifelse(nchar(df$Ref)==1 & nchar(df$Alt)==1, "SNP", "INDEL"),
    Tumor_Sample_Barcode = as.character(sample),
    stringsAsFactors = FALSE
  )
  
  # Map classifications
  maf$Variant_Classification <- case_when(
    maf$Variant_Classification %in% c("nonsynonymous SNV","nonsynonymous_SNV") ~ "Missense_Mutation",
    maf$Variant_Classification == "stopgain" ~ "Nonsense_Mutation",
    maf$Variant_Classification == "stoploss" ~ "Nonstop_Mutation",
    maf$Variant_Classification %in% c("frameshift deletion","frameshift_deletion") ~ "Frame_Shift_Del",
    maf$Variant_Classification %in% c("frameshift insertion","frameshift_insertion") ~ "Frame_Shift_Ins",
    maf$Variant_Classification %in% c("nonframeshift deletion","nonframeshift_deletion") ~ "In_Frame_Del",
    maf$Variant_Classification %in% c("nonframeshift insertion","nonframeshift_insertion") ~ "In_Frame_Ins",
    maf$Variant_Classification == "splicing" ~ "Splice_Site",
    TRUE ~ NA_character_
  )
  
  maf <- maf %>%
    filter(!is.na(Hugo_Symbol),
           !is.na(Variant_Classification))
  
  return(maf)
}

############################################################
# 4. LOAD SOUTH ASIAN MAF
############################################################
files <- list.files("annovar_results", pattern="\\.hg38_multianno\\.txt$", full.names=TRUE)

maf_list <- lapply(files, convert_annovar_to_df)

# 🔥 SAFE BIND (now works)
maf_sa_df <- bind_rows(maf_list)

# Debug genes
cat("\n=== SA RAW GENES ===\n")
print(head(sort(table(maf_sa_df$Hugo_Symbol), decreasing=TRUE), 20))

# Create MAF object
maf_sa <- read.maf(maf = maf_sa_df, vc_nonSyn = NULL)

# Remove noise
maf_sa <- subsetMaf(
  maf = maf_sa,
  query = paste0("!Hugo_Symbol %in% c('", paste(noise_genes, collapse="','"), "')")
)

############################################################
# 5. LOAD TCGA MAF (FIXED)
############################################################

library(data.table)

tcga <- tcgaLoad("HNSC", source = "MC3")

oral_sites <- c("Oral Tongue","Base of tongue","Floor of mouth",
                "Buccal Mucosa","Alveolar Ridge","Hard Palate",
                "Oral Cavity","Lip")

# Filter clinical
eur <- tcga@clinical.data %>%
  dplyr::filter(
    tolower(race) == "white",
    anatomic_neoplasm_subdivision %in% oral_sites
  )

# ✅ Create tcga_eur PROPERLY
tcga_eur <- subsetMaf(
  maf = tcga,
  tsb = eur$Tumor_Sample_Barcode
)

# ✅ Ensure data.table (CRITICAL)
tcga_eur@data <- as.data.table(tcga_eur@data)

# ✅ Remove duplicate columns safely
tcga_eur@data <- tcga_eur@data[, !duplicated(colnames(tcga_eur@data)), with = FALSE]

# ✅ Clean + filter using data.table (NOT dplyr)
tcga_eur@data[, Hugo_Symbol := toupper(as.character(Hugo_Symbol))]

tcga_eur@data <- tcga_eur@data[
  Variant_Classification %in% c(
    "Missense_Mutation","Nonsense_Mutation",
    "Frame_Shift_Del","Frame_Shift_Ins",
    "In_Frame_Del","In_Frame_Ins",
    "Splice_Site","Nonstop_Mutation"
  )
]

# ✅ Apply SAME noise filter
tcga_eur <- subsetMaf(
  maf = tcga_eur,
  query = paste0("!Hugo_Symbol %in% c('", paste(noise_genes, collapse="','"), "')")
)

cat("\n=== SA TOP GENES ===\n")
print(getGeneSummary(maf_sa)[1:15, ])

cat("\n=== TCGA TOP GENES ===\n")
print(getGeneSummary(tcga_eur)[1:15, ])


############################################################
# 🔥 FINAL RECURRENCE FILTER (REMOVE LOW-FREQUENCY NOISE)
############################################################

# Get gene summaries
gs_sa  <- getGeneSummary(maf_sa_filt)
gs_eur <- getGeneSummary(tcga_eur)

# Keep genes mutated in >= 5 samples (adjust if needed)
keep_genes <- union(
  gs_sa$Hugo_Symbol[gs_sa$MutatedSamples >= 5],
  gs_eur$Hugo_Symbol[gs_eur$MutatedSamples >= 5]
)

# Apply to BOTH cohorts
maf_sa_filt <- subsetMaf(
  maf = maf_sa_filt,
  genes = keep_genes
)

tcga_eur <- subsetMaf(
  maf = tcga_eur,
  genes = keep_genes
)

cat("\n=== AFTER RECURRENCE FILTER ===\n")
print(head(getGeneSummary(maf_sa_filt), 10))
print(head(getGeneSummary(tcga_eur), 10))

############################################################
# 🔥 FINAL DRIVER FILTER (SA + TCGA)
############################################################

# Driver gene sets (ensure this is already defined)
# driver_genes <- toupper(c(
#   "TP53","CDKN2A","NOTCH1","FAT1",
#   "CTNNB1","APC","AXIN1","AXIN2","TCF7L2",
#   "PIK3CA","PIK3R1","AKT1","PTEN","MTOR",
#   "KRAS","NRAS","HRAS","BRAF"
# ))

# --------------------------------------------
# Filter South Asian MAF to driver genes only
maf_sa <- subsetMaf(
  maf   = maf_sa,
  genes = driver_genes
)

# Filter TCGA MAF to same driver genes
tcga_eur <- subsetMaf(
  maf   = tcga_eur,
  genes = driver_genes
)

# --------------------------------------------
# Debug / check top genes
cat("\n=== FINAL SA TOP GENES (DRIVERS ONLY) ===\n")
print(getGeneSummary(maf_sa)[1:15, ])

cat("\n=== FINAL TCGA TOP GENES (DRIVERS ONLY) ===\n")
print(getGeneSummary(tcga_eur)[1:15, ])

############################################################
# 5. TMB
############################################################

tmb_sa  <- tmb(maf_sa_filt)
tmb_eur <- tmb(tcga_eur_filt)

tmb_sa$Cohort  <- "South Asian"
tmb_eur$Cohort <- "European"

tmb_df <- dplyr::bind_rows(tmb_sa, tmb_eur)

############################################################
# 6. ROBUST FISHER + EFFECT SIZE
############################################################

get_fisher <- function(m1, m2) {
  
  s1 <- unique(m1@data$Tumor_Sample_Barcode)
  s2 <- unique(m2@data$Tumor_Sample_Barcode)
  
  genes <- union(
    unique(m1@data$Hugo_Symbol),
    unique(m2@data$Hugo_Symbol)
  )
  
  res <- lapply(genes, function(g) {
    
    a  <- length(unique(m1@data$Tumor_Sample_Barcode[m1@data$Hugo_Symbol == g]))
    c2 <- length(unique(m2@data$Tumor_Sample_Barcode[m2@data$Hugo_Symbol == g]))
    
    b <- length(s1) - a
    d <- length(s2) - c2
    
    mat <- matrix(c(a, b, c2, d), 2) + 0.5
    ft  <- fisher.test(mat)
    
    data.frame(
      Gene    = g,
      OR      = as.numeric(ft$estimate),
      CI_low  = ft$conf.int[1],
      CI_high = ft$conf.int[2],
      p       = ft$p.value
    )
  })
  
  df <- dplyr::bind_rows(res)
  df$p_adj <- p.adjust(df$p, "BH")
  df
}

fisher <- get_fisher(maf_sa_filt, tcga_eur_filt)

############################################################
# 7. VOLCANO (PAPER-READY)
############################################################

fisher <- fisher %>%
  dplyr::mutate(
    logOR  = log2(OR),
    neglog = -log10(p_adj)
  )

p_volcano <- ggplot(fisher,
                    aes(logOR, neglog)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(title = "Differential Mutation Enrichment",
       x = "log2(OR)",
       y = "-log10(FDR)")

ggsave(file.path(output_dir,"volcano.pdf"), p_volcano)

############################################################
# 8. FOREST (TOP GENES)
############################################################

driver_genes <- c(
  "TP53","FAT1","NOTCH1","PIK3CA","HRAS","KRAS","CASP8","CDKN2A",
  "NSD1","KMT2D","AJUBA","NFE2L2","KEAP1","PTEN","PIK3R1",
  "CTNNB1","APC","BRAF","NRAS","HRAS","CASP8","TSC1","TSC2",
  "MTOR","EP300","CREBBP","SMAD4","FGFR1","FGFR3"
)

forest <- fisher %>%
  dplyr::filter(
    p_adj < 0.05,
    Gene %in% driver_genes
  ) %>%
  dplyr::arrange(p_adj) %>%
  head(20)

if (nrow(forest) == 0) {
  warning("No driver genes with p_adj < 0.05 for forest plot.")
} else {
  forest$Gene <- factor(forest$Gene, levels = rev(forest$Gene))
  
  p_forest <- ggplot(forest,
                     aes(OR, Gene)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    theme_bw() +
    labs(x = "Odds ratio (log10 scale)",
         y = "Gene",
         title = "Top differential mutation genes (drivers only)")
  
  ggsave(file.path(output_dir,"forest.pdf"), p_forest)
}

############################################################
# 9. KEGG ENRICHMENT
############################################################

top_genes <- fisher %>%
  dplyr::arrange(p_adj) %>%
  head(200) %>%
  dplyr::pull(Gene)

entrez <- bitr(top_genes,
               fromType = "SYMBOL",
               toType   = "ENTREZID",
               OrgDb    = org.Hs.eg.db)

if (!is.null(entrez) && nrow(entrez) > 0) {
  kegg <- enrichKEGG(entrez$ENTREZID, "hsa")
  
  write.csv(as.data.frame(kegg),
            file.path(output_dir,"kegg.csv"),
            row.names = FALSE)
} else {
  warning("No genes mapped to Entrez IDs for KEGG enrichment.")
}

############################################################
# 10. dNdScv (Driver Detection)
############################################################

cat("\n===== Running dNdScv (Driver Detection) =====\n")

if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
remotes::install_github("im3sanger/dndscv")

library(dndscv)

prepare_dndscv <- function(maf_obj) {
  df <- maf_obj@data
  
  dnd <- data.frame(
    sampleID = df$Tumor_Sample_Barcode,
    chr      = gsub("chr", "", df$Chromosome),
    pos      = df$Start_Position,
    ref      = df$Reference_Allele,
    mut      = df$Tumor_Seq_Allele2
  )
  
  dnd <- dnd[complete.cases(dnd), ]
  dnd
}

dnd_sa  <- prepare_dndscv(maf_sa_filt)
dnd_eur <- prepare_dndscv(tcga_eur)

dnd_out_sa  <- dndscv(dnd_sa,  refdb = "hg38")
dnd_out_eur <- dndscv(dnd_eur, refdb = "hg38")

drivers_sa <- dnd_out_sa$sel_cv %>%
  dplyr::filter(qglobal_cv < 0.1)

drivers_eur <- dnd_out_eur$sel_cv %>%
  dplyr::filter(qglobal_cv < 0.1)

write.csv(drivers_sa,
          file.path(output_dir, "dNdScv_SA.csv"),
          row.names = FALSE)

write.csv(drivers_eur,
          file.path(output_dir, "dNdScv_EUR.csv"),
          row.names = FALSE)

cat("dNdScv completed.\n")

############################################################
# 11. MUTATIONAL SIGNATURE ANALYSIS (COSMIC SBS)
############################################################

############################################################
# 11. MUTATIONAL SIGNATURE ANALYSIS (COSMIC SBS, hg38)
############################################################

cat("\n===== Mutational Signature Analysis =====\n")

# 1. Load hg38 genome
library(BSgenome.Hsapiens.UCSC.hg38,
        lib.loc="/media/maazah/Expansion/R_libs/4.5")

# 2. Filter to canonical chromosomes
canonical_chroms <- paste0("chr", c(1:22, "X", "Y"))

maf_sa_clean  <- maf_sa_filt
maf_sa_clean@data  <- maf_sa_clean@data[
  maf_sa_clean@data$Chromosome %in% canonical_chroms, ]

maf_eur_clean <- tcga_eur_filt
maf_eur_clean@data <- maf_eur_clean@data[
  maf_eur_clean@data$Chromosome %in% canonical_chroms, ]

# 3. Recompute trinucleotide matrices (hg38)
tnm_sa <- trinucleotideMatrix(
  maf        = maf_sa_clean,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
)

tnm_eur <- trinucleotideMatrix(
  maf        = maf_eur_clean,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
)

# 4. Extract COSMIC v3 SBS signatures
sig_sa <- extractSignatures(
  mat = tnm_sa$nmf_matrix,
  n   = 2
)

sig_eur <- extractSignatures(
  mat = tnm_eur$nmf_matrix,
  n   = 2
)

# 5. Compare to COSMIC v3 SBS
cosmic_sa  <- compareSignatures(sig_sa$signatures,  sig_db = "SBS")
cosmic_eur <- compareSignatures(sig_eur$signatures, sig_db = "SBS")

# 6. Plot outputs
pdf(file.path(output_dir, "Signatures_SA.pdf"))
plotSignatures(sig_sa)
dev.off()

pdf(file.path(output_dir, "Signatures_EUR.pdf"))
plotSignatures(sig_eur)
dev.off()

pdf(file.path(output_dir, "COSMIC_Similarity_SA.pdf"))
plot(cosmic_sa)
dev.off()

pdf(file.path(output_dir, "COSMIC_Similarity_EUR.pdf"))
plot(cosmic_eur)
dev.off()

cat("Signature analysis completed.\n")


############################################################
# DONE
############################################################

cat("\n===== PUBLICATION PIPELINE COMPLETE =====\n")
