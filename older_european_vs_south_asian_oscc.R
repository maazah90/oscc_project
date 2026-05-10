# TCGA BIO LINKS Install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages(c("fs"), lib = "~/Rlibs", ask = FALSE, update = FALSE)

BiocManager::install("GOSemSim", lib = "/media/maazah/Expansion/R_libs", ask = FALSE, update = FALSE)
BiocManager::install("maftools", lib = "~/Rlibs")
BiocManager::install("Rhtslib")


# ============================================
# 📦 LOAD LIBRARIES
# ============================================
library(TCGAbiolinks)   # for querying TCGA data
library(maftools)       # for mutation data handling
library(dplyr)
library(stringr)
library(ggplot2)

# ============================================
# 📁 PATHS
# ============================================
output_dir <- "maf_results"
dir.create(output_dir, showWarnings = FALSE)

# ============================================
# 🧪 LOAD SOUTH ASIAN MAF (your filtered cohort)
# ============================================
# Assuming maf_filtered exists in environment from prior script
# If not reload it:
# load("maf_filtered.RData")

# Quick check
print(head(getGeneSummary(maf_filtered), 10))

# ============================================
# 🔄 DOWNLOAD TCGA HNSC MAF (MC3 dataset)
# ============================================
# tcgaLoad downloads a curated MAF for the cohort
tcga_hnsc <- tcgaLoad(study = "HNSC", source = "MC3")

# ============================================
# 📌 DOWNLOAD CLINICAL DATA
# ============================================
clin <- GDCquery_clinic("TCGA-HNSC", type = "clinical")

# ============================================
# 📌 CLEAN CLINICAL TABLE
# ============================================
clin <- clin %>%
  mutate(
    Tumor_Sample_Barcode = toupper(substr(submitter_id, 1, 12))
  ) %>%
  select(Tumor_Sample_Barcode, race)

# ============================================
# 🧬 FILTER EUROPEAN ANCESTRY (RACE == WHITE)
# ============================================
white_patients <- clin$Tumor_Sample_Barcode[clin$race == "WHITE"]

# Extract WHITE samples directly from the MAF clinical data
white_samples <- tcga_hnsc@clinical.data %>%
  filter(race == "WHITE") %>%
  pull(Tumor_Sample_Barcode)

# Subset the MAF
tcga_hnsc_eur_oscc <- subsetMaf(
  maf = tcga_hnsc,
  tsb = white_samples
)


# ============================================
# 📊 COMPARISON PLOTS
# ============================================

# -----------------------------
# 1 European TCGA Summary
# -----------------------------
euro_summary <- getGeneSummary(tcga_hnsc_eur_oscc)

write.table(euro_summary,
            file.path(output_dir, "TCGA_European_HNSC_GeneSummary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 2 Side‑by‑side Oncoplots
# -----------------------------
# Top genes in each cohort
top20_south <- getGeneSummary(maf_filtered) %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

top20_euro <- euro_summary %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

# Union of top genes
comparison_genes <- unique(c(top20_south, top20_euro))

# South Asian
pdf(file.path(output_dir, "Comparison_Oncoplot_South_Asian.pdf"), width = 10, height = 8)
oncoplot(maf_filtered, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

# European
pdf(file.path(output_dir, "Comparison_Oncoplot_European_TCGAHNSC.pdf"), width = 10, height = 8)
oncoplot(tcga_hnsc_eur_oscc, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

# -----------------------------
# 3️⃣ Gene Frequency Barplots
# -----------------------------
# Combine summaries
euro_summary$Cohort <- "European"
south_summary <- getGeneSummary(maf_filtered) %>%
  mutate(Cohort = "South Asian")

combined_summary <- bind_rows(
  euro_summary %>% filter(Hugo_Symbol %in% comparison_genes),
  south_summary %>% filter(Hugo_Symbol %in% comparison_genes)
)

# Draw plot
p <- ggplot(combined_summary,
            aes(x = reorder(Hugo_Symbol, -MutatedSamples),
                y = MutatedSamples,
                fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic(base_size = 14) +
  labs(title = "Mutation Frequencies: European vs South Asian",
       x = "",
       y = "Mutated Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Comparison_Mutation_Frequencies.pdf"), p, width = 10, height = 6)

# -----------------------------
# 4️⃣ VAF Distribution Comparison (FIXED)
# -----------------------------

# ============================================
# 🧬 SOUTH ASIAN VAF (USE ORIGINAL TABLE)
# ============================================

vaf1 <- maf_vaf %>%
  filter(
    Variant_Classification %in% c(
      "Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
      "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site"
    )
  ) %>%
  pull(VAF)

# Clean values
vaf1 <- vaf1[!is.na(vaf1) & is.finite(vaf1) & vaf1 > 0 & vaf1 <= 1]

cat("South Asian VAF count:", length(vaf1), "\n")


# ============================================
# 🧬 EUROPEAN TCGA VAF (SAFE)
# ============================================

tcga_data <- tcga_hnsc_eur@data

# 🔥 FIX duplicate column names (CRITICAL)
colnames(tcga_data) <- make.unique(colnames(tcga_data))

# Check allele count columns exist
if(all(c("t_alt_count", "t_ref_count") %in% colnames(tcga_data))){
  
  tcga_data <- tcga_data %>%
    mutate(
      VAF = ifelse(
        (t_alt_count + t_ref_count) > 0,
        t_alt_count / (t_alt_count + t_ref_count),
        NA
      )
    )
  
  vaf2 <- tcga_data$VAF
  
} else {
  stop("❌ TCGA MAF missing allele count columns")
}

# Clean values
vaf2 <- vaf2[!is.na(vaf2) & is.finite(vaf2) & vaf2 > 0 & vaf2 <= 1]

cat("European VAF count:", length(vaf2), "\n")


# ============================================
# 📊 COMBINE DATA
# ============================================

df <- data.frame(
  VAF = c(vaf1, vaf2),
  Cohort = c(
    rep("South Asian", length(vaf1)),
    rep("European", length(vaf2))
  )
)

# Debug check
print(table(df$Cohort))


# ============================================
# 📊 DENSITY PLOT
# ============================================

p_density <- ggplot(df, aes(x = VAF, fill = Cohort)) +
  geom_density(alpha = 0.4) +
  theme_classic(base_size = 14) +
  labs(
    title = "Variant Allele Frequency Distribution",
    x = "Variant Allele Frequency (VAF)",
    y = "Density",
    fill = "Cohort"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "Comparison_VAF_Distribution.pdf"),
  p_density,
  width = 10,
  height = 6
)


# ============================================
# 📊 BOXPLOT (VERY IMPORTANT FOR PRESENTATION)
# ============================================

p_box <- ggplot(df, aes(x = Cohort, y = VAF, fill = Cohort)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.color = "red") +
  theme_classic(base_size = 14) +
  labs(
    title = "VAF Comparison: South Asian vs European Cohorts",
    x = "Cohort",
    y = "Variant Allele Frequency (VAF)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(output_dir, "Comparison_VAF_Boxplot.pdf"),
  p_box,
  width = 8,
  height = 6
)

# -----------------------------
# 5️⃣ Cohort Comparison (BEST)
# -----------------------------
comp <- mafCompare(
  m1 = maf_filtered,
  m2 = tcga_hnsc_eur,
  m1Name = "South Asian",
  m2Name = "European"
)

# Forest plot (VERY GOOD FIGURE)
pdf(file.path(output_dir, "Comparison_MAF_Compare_Forest.pdf"),
    width = 8, height = 10)

forestPlot(comp)

dev.off()
# ============================================
# 🎯 EXPORT ANCESTRY COUNTS
# ============================================
table(tcga_hnsc_eur@clinical.data$race)
write.table(tcga_hnsc_eur@clinical.data,
            file.path(output_dir, "TCGA_European_Clinical.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ Comparison analysis complete!\n")

# Figures

vaf_plot_data <- maf_vaf %>%
  filter(!is.na(VAF)) %>%
  filter(Variant_Classification %in% c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins"
  ))

p <- ggplot(vaf_plot_data,
            aes(x = Variant_Classification, y = VAF, fill = Variant_Classification)) +
  
  geom_violin(trim = FALSE, alpha = 0.6) +
  
  geom_boxplot(width = 0.15, fill = "white") +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "VAF Distribution by Mutation Type",
    x = "Mutation Type",
    y = "Variant Allele Frequency"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(output_dir, "Fig_VAF_by_MutationType.pdf"), p, width = 8, height = 6)















debug_fisher <- function(maf1, maf2) {
  
  g1 <- getGeneSummary(maf1)
  g2 <- getGeneSummary(maf2)
  
  total1 <- as.numeric(maf1@summary[maf1@summary$ID == "Samples", ]$summary)
  total2 <- as.numeric(maf2@summary[maf2@summary$ID == "Samples", ]$summary)
  
  all_genes <- union(g1$Hugo_Symbol, g2$Hugo_Symbol)
  
  for (g in all_genes) {
    
    a <- g1[g1$Hugo_Symbol == g, "MutatedSamples"]
    c <- g2[g2$Hugo_Symbol == g, "MutatedSamples"]
    
    # Force numeric or 0
    a <- if (length(a) == 0) 0 else as.numeric(a[[1]])
    c <- if (length(c) == 0) 0 else as.numeric(c[[1]])
    
    b <- total1 - a
    d <- total2 - c
    
    mat <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)
    
    cat("\n============================\n")
    cat("Gene:", g, "\n")
    print(mat)
    cat("Row sums:", rowSums(mat), " Col sums:", colSums(mat), "\n")
    
    # Try Fisher
    tryCatch({
      fisher.test(mat)
    }, error = function(e) {
      cat("\n❌ ERROR for gene:", g, "\n")
      print(e)
      stop("STOPPING — this is the failing gene.")
    })
  }
}
