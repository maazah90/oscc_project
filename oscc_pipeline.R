# ============================================
# 📦 LIBRARIES
# ============================================
library(data.table)
library(dplyr)
library(stringr)
library(maftools)
library(ggplot2)
library(pheatmap)
library(tidyr)

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

driver_genes <- unique(toupper(c(
  "TP53","CDKN2A","NOTCH1","FAT1",
  wnt_genes, pi3k_genes, ras_genes
)))

# ============================================
# 🔄 ANNOVAR → MAF FUNCTION (SAFE)
# ============================================
convert_annovar_to_maf <- function(file) {
  
  df <- fread(file, data.table = FALSE)
  
  # Ensure required columns exist
  required_cols <- c("Chr","Start","End","Ref","Alt","Gene.ensGene","ExonicFunc.ensGene")
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0){
    stop(paste("Missing columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Fix types
  df$Chr <- as.character(df$Chr)
  df$Start <- as.numeric(df$Start)
  df$End <- as.numeric(df$End)
  df$Ref <- as.character(df$Ref)
  df$Alt <- as.character(df$Alt)
  
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
  
  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_name <- str_replace(sample_name, ".multianno", "")
  
  maf <- data.frame(
    Hugo_Symbol = toupper(str_trim(str_split(df$Gene.ensGene, "[,;]", simplify = TRUE)[,1])),
    Chromosome = df$Chr,
    Start_Position = df$Start,
    End_Position = df$End,
    Reference_Allele = toupper(df$Ref),
    Tumor_Seq_Allele2 = toupper(df$Alt),
    Variant_Classification = df$ExonicFunc.ensGene,
    Variant_Type = case_when(
      nchar(df$Ref) == 1 & nchar(df$Alt) == 1 ~ "SNP",
      nchar(df$Ref) > nchar(df$Alt) ~ "DEL",
      nchar(df$Ref) < nchar(df$Alt) ~ "INS",
      TRUE ~ "SNP"
    ),
    Tumor_Sample_Barcode = sample_name,
    HGVSp_Short = sapply(strsplit(df$AAChange.ensGene, ":"), function(x) tail(x, 1)),
    ClinVar_Pathogenic = df$ClinVar_Pathogenic,
    ClinVar_Significance = df$ClinVar_Significance,
    stringsAsFactors = FALSE
  )
  
  # Standardize variant types
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
  
  return(maf)
}

# ============================================
# 🔗 MERGE
# ============================================
maf_list <- lapply(files, convert_annovar_to_maf)
merged_maf <- bind_rows(maf_list)

# ============================================
# 🧬 ADD VAF FROM AD TABLES
# ============================================
ad_files <- list.files("results/ad_tables", full.names = TRUE)

ad_list <- lapply(ad_files, function(f) {
  
  df <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("Chromosome", "Start_Position", "Ref", "Alt", "AD")
  
  # 🔥 FORCE CONSISTENT TYPES
  df$Chromosome <- as.character(df$Chromosome)
  df$Start_Position <- as.numeric(df$Start_Position)
  df$Ref <- as.character(df$Ref)
  df$Alt <- as.character(df$Alt)
  
  sample <- gsub("_AD.txt", "", basename(f))
  df$Tumor_Sample_Barcode <- sample
  
  return(df)
})

ad_data <- bind_rows(ad_list)

ad_data <- ad_data %>%
  separate(AD, into = c("Ref_Count", "Alt_Count"), sep = ",", convert = TRUE)

# Fix formats
# --------------------------------
# 🔥 FIX ALL MERGE MISMATCH ISSUES
# --------------------------------

# Standardize chromosome format
merged_maf$Chromosome <- gsub("^chr", "", merged_maf$Chromosome, ignore.case = TRUE)
ad_data$Chromosome <- gsub("^chr", "", ad_data$Chromosome, ignore.case = TRUE)

# Ensure same type
merged_maf$Chromosome <- as.character(merged_maf$Chromosome)
ad_data$Chromosome <- as.character(ad_data$Chromosome)

merged_maf$Start_Position <- as.numeric(merged_maf$Start_Position)
ad_data$Start_Position <- as.numeric(ad_data$Start_Position)

# Standardize alleles
merged_maf$Reference_Allele <- toupper(trimws(merged_maf$Reference_Allele))
merged_maf$Tumor_Seq_Allele2 <- toupper(trimws(merged_maf$Tumor_Seq_Allele2))

ad_data$Ref <- toupper(trimws(ad_data$Ref))
ad_data$Alt <- toupper(trimws(ad_data$Alt))

# Standardize sample names (VERY IMPORTANT)
merged_maf$Tumor_Sample_Barcode <- trimws(merged_maf$Tumor_Sample_Barcode)
ad_data$Tumor_Sample_Barcode <- trimws(ad_data$Tumor_Sample_Barcode)

# Fix sample names in merged_maf
merged_maf$Tumor_Sample_Barcode <- gsub("\\.hg38$", "", merged_maf$Tumor_Sample_Barcode)

# Merge
# --------------------------------
# 🔗 ROBUST MERGE
# --------------------------------

maf_vaf <- merge(
  merged_maf,
  ad_data,
  by = c("Chromosome", "Start_Position", "Tumor_Sample_Barcode"),
  all.x = TRUE
)

# --------------------------------
# 🧪 DEBUG
# --------------------------------
cat("Total variants:", nrow(merged_maf), "\n")
cat("Matched AD:", sum(!is.na(maf_vaf$Ref_Count)), "\n")

# --------------------------------
# 🧬 OPTIONAL allele match flag
# --------------------------------
maf_vaf <- maf_vaf %>%
  mutate(
    allele_match = Reference_Allele == Ref & Tumor_Seq_Allele2 == Alt
  )

# --------------------------------
# 📊 SAFE VAF calculation
# --------------------------------
maf_vaf <- maf_vaf %>%
  mutate(
    VAF = ifelse(
      !is.na(Ref_Count) & !is.na(Alt_Count) & (Ref_Count + Alt_Count) > 0,
      Alt_Count / (Ref_Count + Alt_Count),
      NA
    )
  )

cat("VAF calculated:", sum(!is.na(maf_vaf$VAF)), "\n")
summary(maf_vaf$VAF)

# ============================================
# 🧬 CREATE MAF OBJECT
# ============================================
valid_classes <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
                   "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site")

maf_vaf <- maf_vaf %>%
  filter(Variant_Classification %in% valid_classes)

maf <- read.maf(maf = maf_vaf)

# ============================================
# 🧹 FILTER NOISE
# ============================================
noise_genes <- toupper(c(
  "TTN","MUC16","MUC4","MUC5AC","MUC5B",
  "OBSCN","DST","FLG","DNAH1","DNAH2",
  "DNAH5","DNAH8","RYR1","AHNAK2",
  "PCLO","UBR4","LRP1","HMCN2",
  "EPPK1","MUC12","PLEC","SPTBN5",
  "APOB","SBF2","BLTP1","CSMD2"
))

maf_filtered <- subsetMaf(
  maf = maf,
  query = paste0("!Hugo_Symbol %in% c('", paste(noise_genes, collapse="','"), "')"),
  genes = driver_genes
)

# ============================================
# 📊 PATHWAY SUMMARY
# ============================================
pathway_summary <- data.frame(
  Pathway = c("WNT","PI3K_AKT","RAS"),
  Mutated_Samples = c(
    length(unique(maf_filtered@data$Tumor_Sample_Barcode[maf_filtered@data$Hugo_Symbol %in% wnt_genes])),
    length(unique(maf_filtered@data$Tumor_Sample_Barcode[maf_filtered@data$Hugo_Symbol %in% pi3k_genes])),
    length(unique(maf_filtered@data$Tumor_Sample_Barcode[maf_filtered@data$Hugo_Symbol %in% ras_genes]))
  )
)

print(pathway_summary)
# Plot pathway mutation burden
p_pathway <- ggplot(pathway_summary, aes(x = Pathway, y = Mutated_Samples)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic(base_size = 14) +
  labs(title = "Pathway Mutation Burden", y = "Mutated Samples")

ggsave(file.path(output_dir, "Fig_Pathway_Burden.pdf"), p_pathway, width = 6, height = 5)

# ============================================
# 📊 STANDARD PLOTS (unchanged)
# ============================================
gene_summary <- getGeneSummary(maf_filtered)

top_genes <- gene_summary %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

pdf(file.path(output_dir, "Fig_Oncoplot.pdf"), width = 10, height = 8)
oncoplot(maf_filtered, genes = top_genes, removeNonMutated = FALSE)
dev.off()

# ============================================
# 🧪 TP53 CHECK
# ============================================
print(
  maf_vaf %>%
    filter(Hugo_Symbol == "TP53") %>%
    count(Variant_Classification)
)

# ===========================================
# Lollipop plots
# ===========================================
genes_for_lollipop <- unique(c(driver_genes, wnt_genes))

for(gene in genes_for_lollipop){
  
  if(gene %in% maf_filtered@data$Hugo_Symbol){
    
    pdf(file.path(output_dir, paste0("Fig_Lollipop_", gene, ".pdf")), 
        width = 10, height = 4)
    
    lollipopPlot(
      maf = maf_filtered,
      gene = gene,
      AACol = "HGVSp_Short",
      showMutationRate = TRUE
    )
    
    dev.off()
  }
}

#=============================================
# WNT-ONLY ONCOPLOT
#=============================================

wnt_present <- wnt_genes[wnt_genes %in% maf_filtered@data$Hugo_Symbol]

pdf(file.path(output_dir, "Fig_WNT_Oncoplot.pdf"), width = 10, height = 6)

oncoplot(
  maf_filtered,
  genes = wnt_present,
  removeNonMutated = FALSE
)

dev.off()

#=============================================
#  Histogram and BOXPLOT
#=============================================

sample_summary <- getSampleSummary(maf_filtered)

p_hist <- ggplot(sample_summary, aes(x = total)) +
  geom_histogram(bins = 15, fill = "steelblue", color = "black", alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of Mutation Burden per Sample",
    x = "Number of Mutations per Sample",
    y = "Number of Samples"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(output_dir, "Fig_MutationBurden_Hist.pdf"), 
       p_hist, width = 7, height = 5)


p_violin <- ggplot(sample_summary, aes(x = "OSCC", y = total, fill = "OSCC")) +
  
  geom_violin(alpha = 0.6, trim = FALSE) +
  
  geom_boxplot(
    width = 0.15,
    fill = "white",
    outlier.color = "red"
  ) +
  
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Mutation Burden Distribution in OSCC",
    x = "Cohort",
    y = "Number of Mutations per Sample",
    fill = "Cohort"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "Fig_MutationBurden_Violin.pdf"),
  p_violin,
  width = 6,
  height = 5
)

#=============================================
#  Variant Type Proportions
#=============================================

variant_counts <- maf_vaf %>%
  count(Variant_Classification)

ggplot(variant_counts, aes(x = "", y = n, fill = Variant_Classification)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Variant Classification Proportions")

ggsave(file.path(output_dir, "Fig_Variant_Pie.pdf"), width = 6, height = 6)


#=============================================
#  Heat Map 
#=============================================

# --------------------------------
# 1️⃣ Build pathway mutation matrix
# --------------------------------
pathway_counts <- maf_filtered@data %>%
  mutate(
    Pathway = case_when(
      Hugo_Symbol %in% wnt_genes ~ "WNT",
      Hugo_Symbol %in% pi3k_genes ~ "PI3K",
      Hugo_Symbol %in% ras_genes ~ "RAS",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Pathway)) %>%
  group_by(Tumor_Sample_Barcode, Pathway) %>%
  summarise(Mutations = n(), .groups = "drop")

# Convert to matrix
heatmap_df <- tidyr::pivot_wider(
  pathway_counts,
  names_from = Pathway,
  values_from = Mutations,
  values_fill = 0
)

# Set rownames
mat <- as.data.frame(heatmap_df)
rownames(mat) <- mat$Tumor_Sample_Barcode
mat$Tumor_Sample_Barcode <- NULL

mat <- as.matrix(mat)

# --------------------------------
# 2️⃣ Optional: log transform (recommended)
# --------------------------------
mat_log <- log1p(mat)

# --------------------------------
# 3️⃣ Plot clustered heatmap
# --------------------------------

library(pheatmap)

pdf(file.path(output_dir, "Fig_Pathway_Heatmap_Clustered_Clean.pdf"),
    width = 8, height = 10)

pheatmap(
  mat_log,
  
  # Clustering
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  # SHOW LABELS
  show_rownames = FALSE,        # too many samples → hide
  show_colnames = TRUE,         # show pathway names
  
  # MAKE PATHWAYS CLEAR
  fontsize_col = 14,
  angle_col = 0,
  
  # CLEAN COLORS
  color = colorRampPalette(c("lightyellow", "orange", "red"))(100),
  
  # REMOVE GRID CLUTTER
  border_color = "black",
  
  # BETTER SCALING
  breaks = seq(0, max(mat_log), length.out = 100),
  
  # Title
  main = "Pathway Mutation Burden Across OSCC Samples"
)

dev.off()

# =============================================
#  Top Genes
#==============================================

top20 <- gene_summary %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20)

ggplot(top20, aes(x = reorder(Hugo_Symbol, MutatedSamples), y = MutatedSamples)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  labs(title = "Top 20 Mutated Genes", x = "", y = "Samples")

ggsave(file.path(output_dir, "Fig_Top20_Genes.pdf"), width = 6, height = 8)



# =============================================
#  WNT-Summary
#==============================================
wnt_summary <- maf_filtered@data %>%
  filter(Hugo_Symbol %in% wnt_genes) %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    Mutations = n(),
    Samples = n_distinct(Tumor_Sample_Barcode)
  ) %>%
  arrange(desc(Mutations))

write.table(
  wnt_summary,
  file = file.path(output_dir, "WNT_Summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =============================================
#  Top Mutated Genes table
#==============================================

top_genes_table <- gene_summary %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20)

write.table(
  top_genes_table,
  file = file.path(output_dir, "Table_Top_Mutated_Genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =============================================
#  Mutation Burden table
#==============================================

sample_table <- getSampleSummary(maf_filtered)

write.table(
  sample_table,
  file = file.path(output_dir, "Table_Mutation_Burden.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =============================================
#  ClinVar table
#==============================================

clinvar_table <- maf_filtered@data %>%
  filter(grepl("Pathogenic", ClinVar_Significance, ignore.case = TRUE)) %>%
  select(
    Hugo_Symbol,
    Tumor_Sample_Barcode,
    Variant_Classification,
    HGVSp_Short,
    ClinVar_Significance
  )

write.table(
  clinvar_table,
  file = file.path(output_dir, "Table_ClinVar_Pathogenic.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)




















# Lollipop (SA)
sa_lollipop <- maf_sa@data %>%
  filter(Hugo_Symbol %in% driver_genes) %>%
  group_by(Hugo_Symbol) %>%
  summarise(Freq = n_distinct(Tumor_Sample_Barcode), .groups = "drop")

p_lollipop <- ggplot(sa_lollipop,
                     aes(reorder(Hugo_Symbol, Freq), Freq)) +
  geom_segment(aes(xend = Hugo_Symbol, y = 0, yend = Freq), color = "grey70") +
  geom_point(size = 4, color = "#D55E00") +
  coord_flip() +
  theme_pub() +
  labs(title = "South Asian Mutation Spectrum")

ggsave(
  filename = file.path(output_dir, "Fig4A_SA_Lollipop.pdf"),
  plot = p_lollipop,
  width = 7,
  height = 6
)

