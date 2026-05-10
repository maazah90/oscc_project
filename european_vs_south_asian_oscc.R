# ============================================================
# Install TCGA Biolinks if missing
# ============================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

# ============================================================
# Load Libraries
# ============================================================
library(TCGAbiolinks)   # TCGA data access
library(maftools)       # Mutation data analysis
library(dplyr)
library(stringr)
library(ggplot2)

# ============================================================
# Output Directory
# ============================================================
output_dir <- "maf_results"
dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# Load South Asian MAF (previously filtered cohort)
# ============================================================
# Assuming maf_filtered exists in the environment
# If not, load from file:
# load("maf_filtered.RData")

# Quick check
print(head(getGeneSummary(maf_filtered), 10))

# ============================================================
# Download TCGA HNSC MAF (MC3 dataset)
# ============================================================
tcga_hnsc <- tcgaLoad(study = "HNSC", source = "MC3")

# ============================================================
# Get TCGA Clinical Data
# ============================================================
clin_query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Clinical",
  file.type = "xml"
)

GDCdownload(clin_query)
clin <- GDCprepare_clinic(clin_query, clinical.info = "patient")

# Clean clinical table
clin <- clin %>%
  mutate(
    Tumor_Sample_Barcode = toupper(substr(bcr_patient_barcode, 1, 12)),
    race = ifelse(is.na(race), "UNKNOWN", race)
  ) %>%
  select(Tumor_Sample_Barcode, race)

# ============================================================
# Merge TCGA Clinical Metadata into MAF Object
# ============================================================
tcga_hnsc@clinical.data <- tcga_hnsc@clinical.data %>%
  mutate(Tumor_Sample_Barcode = toupper(Tumor_Sample_Barcode)) %>%
  left_join(clin, by = "Tumor_Sample_Barcode")

# ============================================================
# Filter European Ancestry (race == WHITE)
# ============================================================
euro_samples <- tcga_hnsc@clinical.data %>%
  filter(race == "WHITE") %>%
  pull(Tumor_Sample_Barcode)

tcga_hnsc_eur <- subsetMaf(
  maf = tcga_hnsc,
  query = paste0(
    "Tumor_Sample_Barcode %in% c('",
    paste(euro_samples, collapse = "','"),
    "')"
  )
)

# ============================================================
# Comparison Plots
# ============================================================

# ---------------------------
# European TCGA Summary
# ---------------------------
euro_summary <- getGeneSummary(tcga_hnsc_eur)

write.table(
  euro_summary,
  file.path(output_dir, "TCGA_European_HNSC_GeneSummary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---------------------------
# Side-by-side Oncoplots
# ---------------------------

# Top genes in each cohort
top20_south <- getGeneSummary(maf_filtered) %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

top20_euro <- euro_summary %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

# Combined top genes
comparison_genes <- unique(c(top20_south, top20_euro))

# South Asian Oncoplot
pdf(file.path(output_dir, "Comparison_Oncoplot_South_Asian.pdf"),
    width = 10, height = 8)
oncoplot(maf_filtered,
         genes = comparison_genes,
         removeNonMutated = FALSE)
dev.off()

# European Oncoplot
pdf(file.path(output_dir, "Comparison_Oncoplot_European_TCGAHNSC.pdf"),
    width = 10, height = 8)
oncoplot(tcga_hnsc_eur,
         genes = comparison_genes,
         removeNonMutated = FALSE)
dev.off()

# ---------------------------
# Gene Frequency Barplots
# ---------------------------

euro_summary$Cohort <- "European"
south_summary <- getGeneSummary(maf_filtered) %>%
  mutate(Cohort = "South Asian")

combined_summary <- bind_rows(
  euro_summary %>% filter(Hugo_Symbol %in% comparison_genes),
  south_summary %>% filter(Hugo_Symbol %in% comparison_genes)
)

p <- ggplot(combined_summary,
            aes(x = reorder(Hugo_Symbol, -MutatedSamples),
                y = MutatedSamples,
                fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic(base_size = 14) +
  labs(
    title = "Mutation Frequencies: European vs South Asian",
    x = "",
    y = "Mutated Samples"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Comparison_Mutation_Frequencies.pdf"),
       p, width = 10, height = 6)

# ---------------------------
# Variant Classification Summary Comparison
# ---------------------------

pdf(file.path(output_dir, "Comparison_Variant_Classification.pdf"),
    width = 10, height = 6)
plotVafComparison(
  maf1 = maf_filtered,
  maf2 = tcga_hnsc_eur,
  cohort1Name = "South Asian",
  cohort2Name = "European"
)
dev.off()

# ---------------------------
# Co-occurrence Comparison
# ---------------------------

pdf(file.path(output_dir, "Comparison_CoMutation_European_vs_South.pdf"),
    width = 10, height = 8)
coMutSummary(
  maf1 = maf_filtered,
  maf2 = tcga_hnsc_eur,
  cohort1Name = "South Asian",
  cohort2Name = "European"
)
dev.off()

# ============================================================
# Export Ancestry Counts for TCGA White-Only Cohort
# ============================================================
table(tcga_hnsc_eur@clinical.data$race)

write.table(
  tcga_hnsc_eur@clinical.data,
  file.path(output_dir, "TCGA_European_Clinical.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Comparison analysis complete.\n")

############################################################
## 6. Mutation Burden (TMB) Analysis
############################################################

# Calculate TMB for each sample
tmb_sa <- tmb(sa_maf)
tmb_tcga <- tmb(tcga_maf)

# Merge TMB table with sample metadata
tmb_sa <- tmb_sa %>%
  inner_join(sa_metadata, by = c("Tumor_Sample_Barcode" = "Sample_ID"))

tmb_tcga <- tmb_tcga %>%
  inner_join(tcga_metadata, by = c("Tumor_Sample_Barcode" = "Sample_ID"))

# Combine for comparison
tmb_combined <- bind_rows(
  tmb_sa %>% mutate(Cohort = "SouthAsian"),
  tmb_tcga %>% mutate(Cohort = "TCGA_White")
)

# Boxplot
ggplot(tmb_combined, aes(x = Cohort, y = total_perMB)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Tumor Mutational Burden (TMB) Comparison",
    x = "Cohort",
    y = "Mutations per MB"
  )

# Wilcoxon rank-sum test
tmb_wilcox <- wilcox.test(
  total_perMB ~ Cohort,
  data = tmb_combined
)

print(tmb_wilcox)

############################################################
## 7. Statistical Testing: Fisher's Exact Test per Gene
############################################################

# Function to compute gene-wise Fisher test
get_fisher <- function(maf1, maf2, label1, label2) {
  g1 <- getGeneSummary(maf1)
  g2 <- getGeneSummary(maf2)
  
  all_genes <- union(g1$Hugo_Symbol, g2$Hugo_Symbol)
  fisher_results <- data.frame()
  
  for (g in all_genes) {
    a <- g1[g1$Hugo_Symbol == g, "MutatedSamples"]
    if (length(a) == 0) a <- 0
    
    c <- g2[g2$Hugo_Symbol == g, "MutatedSamples"]
    if (length(c) == 0) c <- 0
    
    b <- maf1@summary$Tumor_Samples - a
    d <- maf2@summary$Tumor_Samples - c
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat)
    
    fisher_results <- rbind(
      fisher_results,
      data.frame(
        Gene = g,
        OR = ft$estimate,
        p = ft$p.value,
        label1 = label1,
        label2 = label2
      )
    )
  }
  
  fisher_results$p_adj <- p.adjust(fisher_results$p, method = "BH")
  fisher_results
}

fisher_table <- get_fisher(sa_maf, tcga_maf, "SouthAsian", "TCGA_White")

write.csv(fisher_table, "results/fisher_gene_level_results.csv", row.names = FALSE)

############################################################
## 8. UpSet Plots for Gene Intersection
############################################################

upsetGenes(sa_maf, tcga_maf, 
           setNames = c("South Asian", "TCGA White"), 
           textScale = 1.2)

############################################################
## 9. Volcano Plot of Gene Differences (Using Fisher Results)
############################################################

fisher_table$logOR <- log2(fisher_table$OR)
fisher_table$neglog10p <- -log10(fisher_table$p_adj)

ggplot(fisher_table, aes(x = logOR, y = neglog10p)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Volcano Plot – Differential Mutation Frequency",
    x = "log2(Odds Ratio)",
    y = "-log10(adjusted p-value)"
  )

############################################################
## 10. PCA of Binary Mutation Matrix
############################################################

# Convert MAF → wide mutation matrix
binary_matrix <- createOncoMatrix(sa_tcga_combined)$numericMatrix
# NA → 0
binary_matrix[is.na(binary_matrix)] <- 0

# PCA
pca_res <- prcomp(t(binary_matrix), scale. = TRUE)

# Extract PC scores
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

# Annotate cohort
pca_df <- pca_df %>%
  left_join(
    rbind(
      sa_metadata %>% mutate(Cohort = "SouthAsian"),
      tcga_metadata %>% mutate(Cohort = "TCGA_White")
    ),
    by = c("Sample" = "Sample_ID")
  )

# PCA Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCA of Binary Mutation Matrix",
    x = "PC1",
    y = "PC2"
  )

############################################################
## 11. KEGG Pathway Enrichment (ClusterProfiler)
############################################################

# Top mutated genes per cohort
top_sa <- getGeneSummary(sa_maf) %>%
  dplyr::arrange(desc(MutatedSamples)) %>%
  dplyr::slice(1:200) %>%
  pull(Hugo_Symbol)

top_tcga <- getGeneSummary(tcga_maf) %>%
  dplyr::arrange(desc(MutatedSamples)) %>%
  dplyr::slice(1:200) %>%
  pull(Hugo_Symbol)

# Convert to Entrez IDs
gene2entrez <- bitr(
  c(top_sa, top_tcga),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_sa <- gene2entrez %>% filter(SYMBOL %in% top_sa) %>% pull(ENTREZID)
entrez_tcga <- gene2entrez %>% filter(SYMBOL %in% top_tcga) %>% pull(ENTREZID)

# KEGG
kegg_sa <- enrichKEGG(gene = entrez_sa, organism = "hsa")
kegg_tcga <- enrichKEGG(gene = entrez_tcga, organism = "hsa")

# Save tables
write.csv(as.data.frame(kegg_sa), "results/kegg_sa.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_tcga), "results/kegg_tcga.csv", row.names = FALSE)

# Dotplots
dotplot(kegg_sa) + ggtitle("KEGG – South Asian")
dotplot(kegg_tcga) + ggtitle("KEGG – TCGA White")

############################################################
## 12. COSMIC SBS Mutational Signatures
############################################################

# Extract trinucleotide matrices
sa_mat <- trinucleotideMatrix(maf = sa_maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
tcga_mat <- trinucleotideMatrix(maf = tcga_maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

# Compare with COSMIC v3
sa_fit <- fitToSignatures(sa_mat, get_known_signatures())
tcga_fit <- fitToSignatures(tcga_mat, get_known_signatures())

# Plot signature contributions
plotSignatureContribution(sa_fit$contribution, title = "South Asian – SBS Signatures")
plotSignatureContribution(tcga_fit$contribution, title = "TCGA White – SBS Signatures")

############################################################
## END MESSAGE 2
############################################################

############################################################
## 13. Utility: Export All Plots (Optional)
############################################################
# You can uncomment this if you want automatic plot saving

# ggsave("results/tmb_boxplot.png", width = 6, height = 5, dpi = 300)
# ggsave("results/volcano_plot.png", width = 6, height = 5, dpi = 300)
# ggsave("results/pca_plot.png", width = 6, height = 5, dpi = 300)

############################################################
## 14. Save Core R Objects for Reproducibility
############################################################

save(
  sa_maf, tcga_maf,
  sa_metadata, tcga_metadata,
  fisher_table,
  tmb_combined,
  kegg_sa, kegg_tcga,
  sa_mat, tcga_mat,
  sa_fit, tcga_fit,
  file = "results/all_objects.RData"
)

############################################################
## 15. Summary Output to Console
############################################################

cat("\n===== PIPELINE COMPLETED SUCCESSFULLY =====\n\n")

cat("South Asian samples:\n")
print(summary(sa_maf))

cat("\nTCGA White OSCC samples:\n")
print(summary(tcga_maf))

cat("\nTop Significant Genes (Fisher, adj p < 0.05):\n")
print(
  fisher_table %>%
    dplyr::filter(p_adj < 0.05) %>%
    dplyr::arrange(p_adj) %>%
    dplyr::slice(1:20)
)

cat("\nKEGG Top Pathways – South Asian:\n")
print(head(as.data.frame(kegg_sa)[,1:5]))

cat("\nKEGG Top Pathways – TCGA White:\n")
print(head(as.data.frame(kegg_tcga)[,1:5]))

cat("\n===== END OF SCRIPT =====\n")

############################################################
## 16. Session Info for Full Reproducibility
############################################################

sessionInfo()