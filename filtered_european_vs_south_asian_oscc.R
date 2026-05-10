# ============================================
# 📦 LOAD LIBRARIES
# ============================================
library(maftools)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# ============================================
# 📁 PATHS
# ============================================
output_dir <- "maf_results"
dir.create(output_dir, showWarnings = FALSE)

# ============================================
# 🧪 SOUTH ASIAN MAF (WES cohort)
# ============================================
# maf_filtered: filtered MAF object from South Asian samples
# maf_vaf: MAF data frame with VAFs

# --------------------------------------------
# Check top genes
gene_summary_south <- getGeneSummary(maf_filtered)
head(gene_summary_south, 10)

# ============================================
# 🔄 EUROPEAN TCGA HNSC (MC3 MAF + clinical)
# ============================================

tcga_hnsc <- tcgaLoad(study = "HNSC", source = "MC3")
clin_mc3  <- tcga_hnsc@clinical.data

# Exact oral cavity subsites from your table
oral_sites <- c(
  "Oral Tongue",
  "Base of tongue",
  "Floor of mouth",
  "Buccal Mucosa",
  "Alveolar Ridge",
  "Hard Palate",
  "Oral Cavity",
  "Lip"
)

# Optional: check race values
print(table(clin_mc3$race))

euro_oscc_clin <- clin_mc3 %>%
  filter(
    tolower(race) == "white",
    anatomic_neoplasm_subdivision %in% oral_sites
  )

cat("Rows in euro_oscc_clin:", nrow(euro_oscc_clin), "\n")

# Use Tumor_Sample_Barcode directly
matched_samples <- euro_oscc_clin$Tumor_Sample_Barcode

cat("Matched TCGA OSCC samples:", length(matched_samples), "\n")

tcga_hnsc_eur_oscc <- subsetMaf(
  maf = tcga_hnsc,
  tsb = matched_samples
)

cat("Final TCGA European OSCC samples:",
    length(unique(tcga_hnsc_eur_oscc@data$Tumor_Sample_Barcode)), "\n")




# --------------------------------
# Oncoplots
# --------------------------------
pdf(file.path(output_dir, "Comparison_Oncoplot_South_Asian.pdf"), width = 10, height = 8)
oncoplot(maf_filtered, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

pdf(file.path(output_dir, "Comparison_Oncoplot_European_TCGA_OSCC.pdf"), width = 10, height = 8)
oncoplot(tcga_hnsc_eur_oscc, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

# --------------------------------
# Gene Frequency Barplot
# --------------------------------
gene_summary_euro <- getGeneSummary(tcga_hnsc_eur_oscc) %>% mutate(Cohort = "European")
gene_summary_south <- getGeneSummary(maf_filtered) %>% mutate(Cohort = "South Asian")

combined_summary <- bind_rows(
  gene_summary_euro %>% filter(Hugo_Symbol %in% comparison_genes),
  gene_summary_south %>% filter(Hugo_Symbol %in% comparison_genes)
)

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

#############################################################
## VAF Comparison – South Asian vs European
############################################################

library(ggplot2)
library(dplyr)

# Combine VAFs into tidy data frame
df_vaf <- data.frame(
  VAF = c(maf_vaf$VAF, tcga_hnsc_eur_oscc@data$VAF),
  Cohort = c(
    rep("South Asian", length(maf_vaf$VAF)),
    rep("European", length(tcga_hnsc_eur_oscc@data$VAF))
  )
)

# Remove NAs
df_vaf <- df_vaf %>% filter(!is.na(VAF))

# Plot
p_vaf <- ggplot(df_vaf, aes(x = Cohort, y = VAF, fill = Cohort)) +
  geom_violin(alpha = 0.6, trim = FALSE, color = "black") +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.6, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
  scale_fill_manual(values = c("South Asian" = "#E69F00", "European" = "#56B4E9")) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw(base_size = 14) +
  labs(
    title = "Variant Allele Frequency (VAF) Distribution",
    x = "",
    y = "VAF"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none"
  )

# Save
ggsave(file.path(output_dir, "Comparison_VAF_Distribution_Presentation.pdf"),
       p_vaf, width = 7, height = 5)


############################################################
## 6. Tumor Mutational Burden (TMB)
############################################################

# Calculate TMB for each cohort
tmb_sa   <- tmb(maf_filtered)
tmb_tcga <- tmb(tcga_hnsc_eur_oscc)

# Add cohort labels
tmb_sa$Cohort   <- "SouthAsian"
tmb_tcga$Cohort <- "TCGA_White"

# Combine
tmb_combined <- bind_rows(tmb_sa, tmb_tcga)

# Boxplot
ggplot(tmb_combined, aes(x = Cohort, y = total_perMB, fill = Cohort)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Tumor Mutational Burden (TMB) Comparison",
    x = "Cohort",
    y = "Mutations per MB"
  )

# Wilcoxon test
tmb_wilcox <- wilcox.test(total_perMB ~ Cohort, data = tmb_combined)
print(tmb_wilcox)

aggregate(total_perMB ~ Cohort, df, median)


library(ggpubr)

ggboxplot(
  tmb_combined,
  x = "Cohort",
  y = "total_perMB",
  fill = "Cohort",
  palette = "jco"
) +
  stat_compare_means(method = "wilcox.test") +
  theme_bw(base_size = 14) +
  ylab("TMB (mutations per MB)") +
  ggtitle("Tumor Mutation Burden by Cohort")



############################################################
## 12. COSMIC SBS Signatures
############################################################

sa_mat   <- trinucleotideMatrix(maf = maf_filtered, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
tcga_mat <- trinucleotideMatrix(maf = tcga_hnsc_eur_oscc, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

sa_fit   <- fitToSignatures(sa_mat, get_known_signatures())
tcga_fit <- fitToSignatures(tcga_mat, get_known_signatures())

plotSignatureContribution(sa_fit$contribution,   title = "South Asian – SBS Signatures")
plotSignatureContribution(tcga_fit$contribution, title = "TCGA White – SBS Signatures")



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
  fisher_table,
  tmb_combined,
  kegg_sa, kegg_tcga,
  sa_mat, sa_fit,
  tcga_mat, tcga_fit,
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




# ============================================
# 📊 Presentation-Ready Graphs & Tables
# ============================================

# --------------------------------------------
# Updated Fisher function (fixed)
# --------------------------------------------
get_fisher <- function(maf1, maf2, label1, label2) {
  
  g1 <- getGeneSummary(maf1)
  g2 <- getGeneSummary(maf2)
  
  total1 <- as.numeric(maf1@summary$Tumor_Samples)
  total2 <- as.numeric(maf2@summary$Tumor_Samples)
  
  all_genes <- union(g1$Hugo_Symbol, g2$Hugo_Symbol)
  fisher_results <- data.frame()
  
  for (g in all_genes) {
    a <- g1[g1$Hugo_Symbol == g, "MutatedSamples"]
    a <- if (length(a) == 0) 0 else as.numeric(a[[1]])
    
    c <- g2[g2$Hugo_Symbol == g, "MutatedSamples"]
    c <- if (length(c) == 0) 0 else as.numeric(c[[1]])
    
    b <- total1 - a
    d <- total2 - c
    
    vals <- c(a, b, c, d)
    if (any(is.na(vals)) || any(vals < 0) || any(!is.finite(vals))) next
    
    mat <- matrix(vals, nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat)
    
    fisher_results <- rbind(
      fisher_results,
      data.frame(
        Gene = g,
        OR = as.numeric(ft$estimate),
        p = ft$p.value
      )
    )
  }
  
  fisher_results$p_adj <- p.adjust(fisher_results$p, method = "BH")
  fisher_results
}

# ============================================
# 1 Oncoplots
# ============================================
pdf(file.path(output_dir, "Oncoplot_SouthAsian.pdf"), width = 10, height = 8)
oncoplot(maf_filtered, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

pdf(file.path(output_dir, "Oncoplot_European.pdf"), width = 10, height = 8)
oncoplot(tcga_hnsc_eur_oscc, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

# ============================================
# 2 Gene Frequency Barplot
# ============================================
gene_summary_sa <- getGeneSummary(maf_filtered) %>% mutate(Cohort = "South Asian")
gene_summary_eur <- getGeneSummary(tcga_hnsc_eur_oscc) %>% mutate(Cohort = "European")

combined_summary <- bind_rows(
  gene_summary_sa %>% filter(Hugo_Symbol %in% comparison_genes),
  gene_summary_eur %>% filter(Hugo_Symbol %in% comparison_genes)
)

p_freq <- ggplot(combined_summary,
                 aes(x = reorder(Hugo_Symbol, -MutatedSamples),
                     y = MutatedSamples, fill = Cohort)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  theme_classic(base_size = 14) +
  labs(title = "Mutation Frequencies: South Asian vs European",
       x = "", y = "Mutated Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "Mutation_Frequency_Comparison.pdf"), p_freq, width = 10, height = 6)

# ============================================
# 3️⃣ VAF Distribution
# ============================================
df_vaf <- data.frame(
  VAF = c(maf_vaf$VAF, tcga_hnsc_eur_oscc@data$VAF),
  Cohort = c(rep("South Asian", length(maf_vaf$VAF)),
             rep("European", length(tcga_hnsc_eur_oscc@data$VAF)))
)

p_vaf <- ggplot(df_vaf, aes(x = Cohort, y = VAF, fill = Cohort)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.color = "red") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  theme_classic(base_size = 14) +
  labs(title = "Variant Allele Frequency (VAF) Comparison", y = "VAF", x = "") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "VAF_Comparison.pdf"), p_vaf, width = 7, height = 5)

# ============================================
# 4 Scatter Plot: Gene Frequency Comparison
# ============================================
# Clean gene names first
maf_filtered@data$Hugo_Symbol <- toupper(trimws(maf_filtered@data$Hugo_Symbol))
comparison_genes <- toupper(trimws(comparison_genes))

# Robust frequency calculation
calc_freq <- function(maf, genes) {
  total <- length(unique(maf@data$Tumor_Sample_Barcode))
  
  # Compute mutated samples
  freq_df <- maf@data %>%
    filter(Hugo_Symbol %in% genes) %>%
    group_by(Hugo_Symbol) %>%
    summarise(Samples = n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
    mutate(Frequency = Samples / total)
  
  # Ensure all genes in list are present (fill 0 for missing)
  missing_genes <- setdiff(genes, freq_df$Hugo_Symbol)
  if(length(missing_genes) > 0){
    freq_df <- bind_rows(
      freq_df,
      data.frame(Hugo_Symbol = missing_genes,
                 Samples = 0,
                 Frequency = 0)
    )
  }
  
  # Sort descending
  freq_df <- freq_df %>% arrange(desc(Frequency))
  return(freq_df)
}

# Run safely
south_freq <- calc_freq(maf_filtered, comparison_genes) %>% rename(South = Frequency)
euro_freq  <- calc_freq(tcga_hnsc_eur_oscc, comparison_genes) %>% rename(Europe = Frequency)

# Merge for scatter
freq_df <- merge(south_freq, euro_freq, by = "Hugo_Symbol", all = TRUE)
freq_df[is.na(freq_df)] <- 0

freq_df$Pathway <- "Other"
freq_df$Pathway[freq_df$Hugo_Symbol %in% wnt_genes] <- "WNT"
freq_df$Pathway[freq_df$Hugo_Symbol %in% pi3k_genes] <- "PI3K"
freq_df$Pathway[freq_df$Hugo_Symbol %in% ras_genes] <- "RAS"

library(ggrepel)
p_scatter <- ggplot(freq_df, aes(x = Europe, y = South, color = Pathway, label = Hugo_Symbol)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "pink") +
  geom_text_repel(max.overlaps = 20, size = 3.5) +
  scale_color_manual(values = c("WNT" = "#1b9e77", "PI3K" = "#d95f02", "RAS" = "#7570b3", "Other" = "pink")) +
  theme_classic(base_size = 14) +
  labs(title = "Mutation Frequency Comparison: South Asian vs European",
       x = "European OSCC Frequency",
       y = "South Asian OSCC Frequency",
       color = "Pathway") +
  xlim(0,1) + ylim(0,1)

ggsave(file.path(output_dir, "Scatter_Frequency_Comparison.pdf"), p_scatter, width = 8, height = 6)

# ============================================
# 5 Tumor Mutational Burden (TMB)
# ============================================
tmb_sa   <- tmb(maf_filtered)
tmb_eur  <- tmb(tcga_hnsc_eur_oscc)

tmb_sa$Cohort <- "South Asian"
tmb_eur$Cohort <- "European"
tmb_combined <- bind_rows(tmb_sa, tmb_eur)

p_tmb <- ggplot(tmb_combined, aes(x = Cohort, y = total_perMB, fill = Cohort)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "red", outlier.size = 2) +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  theme_bw(base_size = 14) +
  labs(title = "Tumor Mutational Burden (TMB) Comparison", x = "", y = "Mutations per MB") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "TMB_Comparison.pdf"), p_tmb, width = 6, height = 5)


## FISHER PLOT 
get_fisher <- function(maf1, maf2, label1, label2) {
  
  g1 <- getGeneSummary(maf1)
  g2 <- getGeneSummary(maf2)
  
  total1 <- as.numeric(maf1@summary[maf1@summary$ID == "Samples", ]$summary)
  total2 <- as.numeric(maf2@summary[maf2@summary$ID == "Samples", ]$summary)
  
  all_genes <- union(g1$Hugo_Symbol, g2$Hugo_Symbol)
  fisher_table <- data.frame()
  
  for (g in all_genes) {
    
    # Extract mutated counts safely
    a <- g1[g1$Hugo_Symbol == g, "MutatedSamples"]
    a <- if (length(a) == 0) 0 else as.numeric(a)
    
    c <- g2[g2$Hugo_Symbol == g, "MutatedSamples"]
    c <- if (length(c) == 0) 0 else as.numeric(c)
    
    # Compute non-mutated
    b <- total1 - a
    d <- total2 - c
    
    vals <- c(a, b, c, d)
    
    # Skip invalid
    if (any(is.na(vals)) || any(vals < 0) || any(!is.finite(vals))) next
    
    mat <- matrix(vals, nrow = 2, byrow = TRUE)
    
    # Skip degenerate matrices
    if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) next
    
    ft <- fisher.test(mat)
    
    fisher_table <- rbind(
      fisher_table,
      data.frame(
        Gene = g,
        OR   = as.numeric(ft$estimate),
        p    = ft$p.value
      )
    )
  }
  
  fisher_table$p_adj <- p.adjust(fisher_table$p, method = "BH")
  return(fisher_table)
}



# Run Fisher test
fisher_table <- get_fisher(
  maf1 = maf_filtered,
  maf2 = tcga_hnsc_eur_oscc,
  label1 = "SouthAsian",
  label2 = "TCGA_White"
)

# Remove OR <= 0 or non-finite
fisher_table <- fisher_table[
  is.finite(fisher_table$OR) & fisher_table$OR > 0,
]

# Add log-transformed values
fisher_table$logOR      <- log2(fisher_table$OR)
fisher_table$neglog10p  <- -log10(fisher_table$p_adj)

# Save results
write.csv(fisher_table,
          file.path(output_dir, "Fisher_Test_Results.csv"),
          row.names = FALSE)

# Volcano plot
p_volcano <- ggplot(fisher_table, aes(x = logOR, y = neglog10p, color = p_adj < 0.05)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot – Differential Mutation Frequency",
    x = "log2(Odds Ratio)",
    y = "-log10(adjusted p-value)",
    color = "Significant"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "Volcano_Fisher_Improved.pdf"),
       p_volcano, width = 7, height = 5)


# ---------------------------
# 3. Variant Classification Pie Charts
# ---------------------------
plot_variant_pie <- function(maf_obj, cohort_name) {
  var_summary <- maf_obj@data %>% count(Variant_Classification) %>% mutate(Percent = n / sum(n) * 100)
  p <- ggplot(var_summary, aes(x = "", y = Percent, fill = Variant_Classification)) +
    geom_col(color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = paste0(round(Percent, 1), "%")), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = paste("Variant Classification –", cohort_name), fill = "Mutation Type") +
    theme_void(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
  return(p)
}
save_plot(plot_variant_pie(maf_filtered, "South Asian"), "VariantClassification_SouthAsian.pdf")
save_plot(plot_variant_pie(tcga_hnsc_eur_oscc, "TCGA White"), "VariantClassification_TCGAWhite.pdf")


############################################################
## PCA – Improved for Presentation
############################################################

library(ggrepel)  # For non-overlapping labels
library(RColorBrewer)

# Compute PCA
if (ncol(binary_matrix) >= 2) {
  pca_res <- prcomp(t(binary_matrix), scale. = TRUE)
  
  # Variance explained
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  
  pca_df <- data.frame(
    Sample = rownames(pca_res$x),
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Cohort = ifelse(rownames(pca_res$x) %in% maf_filtered@data$Tumor_Sample_Barcode,
                    "South Asian", "European")
  )
  
  # Plot with nicer colors and shapes
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cohort, shape = Cohort)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 14) +
    labs(
      title = "PCA of Binary Mutation Matrix",
      x = paste0("PC1 (", round(var_explained[1]*100,1), "%)"),
      y = paste0("PC2 (", round(var_explained[2]*100,1), "%)")
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "top"
    ) +
    geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 5)
  
  # Save figure
  ggsave(file.path(output_dir, "PCA_BinaryMutation.pdf"), p_pca, width = 8, height = 6)
  
  print(p_pca)
}

# ============================================
# 8. Lollipop plots for key pathways
# ============================================
genes_for_lollipop <- unique(c(wnt_genes, pi3k_genes, ras_genes, driver_genes))

for(gene in genes_for_lollipop){
  if(gene %in% maf_filtered@data$Hugo_Symbol){
    pdf(file.path(output_dir, paste0("Lollipop_", gene, ".pdf")), width = 10, height = 4)
    lollipopPlot(maf_filtered, gene = gene, AACol = "HGVSp_Short", showMutationRate = TRUE)
    dev.off()
  }
}




# ============================================
# 9. Save Top Genes & Pathway Tables
# ============================================
write.csv(getGeneSummary(maf_filtered), file.path(output_dir, "TopGenes_SouthAsian.csv"), row.names = FALSE)
write.csv(getGeneSummary(tcga_hnsc_eur_oscc), file.path(output_dir, "TopGenes_European.csv"), row.names = FALSE)

# Pathway mutation table
pathway_summary <- data.frame(
  Pathway = c("WNT","PI3K","RAS"),
  Mutated_Samples_SouthAsian = sapply(list(wnt_genes, pi3k_genes, ras_genes),
                                      function(x) length(unique(maf_filtered@data$Tumor_Sample_Barcode[maf_filtered@data$Hugo_Symbol %in% x]))),
  Mutated_Samples_European = sapply(list(wnt_genes, pi3k_genes, ras_genes),
                                    function(x) length(unique(tcga_hnsc_eur_oscc@data$Tumor_Sample_Barcode[tcga_hnsc_eur_oscc@data$Hugo_Symbol %in% x])))
)
write.csv(pathway_summary, file.path(output_dir, "PathwayMutationTable.csv"), row.names = FALSE)



############################################################
## KEGG Pathway Enrichment (Fixed)
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# -----------------------------
# Top genes per cohort
# -----------------------------
top_sa <- getGeneSummary(maf_filtered) %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 200) %>%
  pull(Hugo_Symbol)

top_tcga <- getGeneSummary(tcga_hnsc_eur_oscc) %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 200) %>%
  pull(Hugo_Symbol)

# -----------------------------
# SYMBOL → ENTREZ conversion
# -----------------------------
gene2entrez <- bitr(
  unique(c(top_sa, top_tcga)),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Filter valid ENTREZ IDs
entrez_sa <- gene2entrez %>%
  filter(SYMBOL %in% top_sa) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

entrez_tcga <- gene2entrez %>%
  filter(SYMBOL %in% top_tcga) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

# -----------------------------
# Run KEGG enrichment (safe)
# -----------------------------
kegg_sa <- enrichKEGG(
  gene = entrez_sa,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_tcga <- enrichKEGG(
  gene = entrez_tcga,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# -----------------------------
# Save results
# -----------------------------
write.csv(as.data.frame(kegg_sa),
          file.path(output_dir, "KEGG_SA.csv"), row.names = FALSE)

write.csv(as.data.frame(kegg_tcga),
          file.path(output_dir, "KEGG_TCGA.csv"), row.names = FALSE)

# -----------------------------
# Optional: plot top pathways
# -----------------------------
if(nrow(as.data.frame(kegg_sa)) > 0){
  barplot(kegg_sa, showCategory = 10, title = "South Asian KEGG Pathways")
}
if(nrow(as.data.frame(kegg_tcga)) > 0){
  barplot(kegg_tcga, showCategory = 10, title = "European KEGG Pathways")
}




library(clusterProfiler)
library(ggplot2)
library(dplyr)

# Convert KEGG enrichment result to data frame
kegg_df <- as.data.frame(kegg_sa)

# Optional: filter to canonical cancer signaling pathways for presentation
cancer_keywords <- c("PI3K", "AKT", "WNT", "RAS", "MAPK", "p53", "Cell cycle", "Apoptosis")
kegg_df <- kegg_df %>%
  filter(grepl(paste(cancer_keywords, collapse="|"), Description, ignore.case = TRUE))

# Reorder by significance
kegg_df <- kegg_df %>%
  arrange(p.adjust) %>%
  mutate(Description = factor(Description, levels = Description))

# Make dotplot
p_kegg <- ggplot(kegg_df, aes(x = -log10(p.adjust), y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_bw(base_size = 14) +
  labs(
    title = "KEGG Pathway Enrichment – South Asian OSCC",
    x = "-log10(adjusted p-value)",
    y = "",
    size = "Mutated Genes",
    color = "-log10(adjusted p-value)"
  ) +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(output_dir, "KEGG_Clean_SA_OSCC.pdf"), p_kegg, width = 8, height = 6)




library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ----------------------------
# Convert KEGG results to data frames
# ----------------------------
kegg_sa_df <- as.data.frame(kegg_sa) %>%
  mutate(Cohort = "South Asian")

kegg_tcga_df <- as.data.frame(kegg_tcga) %>%
  mutate(Cohort = "European")

# ----------------------------
# Optional: filter for key cancer pathways
# ----------------------------
cancer_keywords <- c("PI3K", "AKT", "WNT", "RAS", "MAPK", "p53", "Cell cycle", "Apoptosis")

kegg_sa_df <- kegg_sa_df %>% filter(grepl(paste(cancer_keywords, collapse = "|"), Description, ignore.case = TRUE))
kegg_tcga_df <- kegg_tcga_df %>% filter(grepl(paste(cancer_keywords, collapse = "|"), Description, ignore.case = TRUE))

# ----------------------------
# Combine both cohorts
# ----------------------------
kegg_combined <- bind_rows(kegg_sa_df, kegg_tcga_df)

# ----------------------------
# Reorder pathways for plotting
# ----------------------------
kegg_combined <- kegg_combined %>%
  arrange(Cohort, p.adjust) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# ----------------------------
# Side-by-side dotplot
# ----------------------------
p_kegg_comparison <- ggplot(kegg_combined, aes(x = -log10(p.adjust), y = Description, color = Cohort, size = Count)) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("South Asian" = "steelblue", "European" = "firebrick")) +
  theme_bw(base_size = 14) +
  labs(
    title = "KEGG Pathway Enrichment Comparison: South Asian vs European OSCC",
    x = "-log10(adjusted p-value)",
    y = "",
    color = "Cohort",
    size = "Mutated Genes"
  ) +
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



# ============================================================
# VAF Density Plot (Normal Distribution Style)
# ============================================================

# Ensure European VAF exists
if (!"VAF" %in% colnames(tcga_hnsc_eur_oscc@data) ||
    all(is.na(tcga_hnsc_eur_oscc@data$VAF))) {
  
  if (all(c("t_alt_count", "t_ref_count") %in% colnames(tcga_hnsc_eur_oscc@data))) {
    tcga_hnsc_eur_oscc@data$VAF <- with(
      tcga_hnsc_eur_oscc@data,
      ifelse(
        (t_alt_count + t_ref_count) > 0,
        t_alt_count / (t_alt_count + t_ref_count),
        NA
      )
    )
  }
}

# Build tidy VAF dataframe
df_vaf <- data.frame(
  VAF = c(maf_vaf$VAF, tcga_hnsc_eur_oscc@data$VAF),
  Cohort = c(
    rep("South Asian", length(maf_vaf$VAF)),
    rep("European",   length(tcga_hnsc_eur_oscc@data$VAF))
  ),
  stringsAsFactors = FALSE
)

df_vaf <- df_vaf[!is.na(df_vaf$VAF), ]

library(ggplot2)

p_vaf_density <- ggplot(df_vaf, aes(x = VAF, fill = Cohort, color = Cohort)) +
  geom_density(alpha = 0.35, linewidth = 1.1) +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  scale_color_manual(values = c("South Asian" = "#1b9e77", "European" = "#d95f02")) +
  theme_bw(base_size = 15) +
  labs(
    title = "Variant Allele Frequency Distribution",
    x = "Variant Allele Frequency (VAF)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "VAF_Density_Comparison.pdf"),
  p_vaf_density, width = 8, height = 5.5
)

