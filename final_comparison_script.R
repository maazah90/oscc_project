############################################################
# OSCC Comparative Genomics Pipeline
# South Asian (WES) vs European (TCGA HNSC OSCC)
############################################################

## 0. Libraries -------------------------------------------------------------

library(maftools)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)
library(UpSetR)

## 1. Paths & helpers -------------------------------------------------------

output_dir <- "maf_results"
dir.create(output_dir, showWarnings = FALSE)

save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(file.path(output_dir, filename), p, width = width, height = height)
}

## 2. Load / define cohorts -------------------------------------------------
# Assumed already in environment:
# - maf_filtered: South Asian MAF object
# - maf_vaf: data.frame with column VAF
# - comparison_genes: character vector of genes
# - wnt_genes, pi3k_genes, ras_genes, driver_genes: character vectors

## 2.1 TCGA HNSC – filter to White OSCC ------------------------------------

tcga_hnsc <- tcgaLoad(study = "HNSC", source = "MC3")
clin_mc3  <- tcga_hnsc@clinical.data

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

euro_oscc_clin <- clin_mc3 %>%
  filter(
    tolower(race) == "white",
    anatomic_neoplasm_subdivision %in% oral_sites
  )

matched_samples <- euro_oscc_clin$Tumor_Sample_Barcode

tcga_hnsc_eur_oscc <- subsetMaf(
  maf = tcga_hnsc,
  tsb = matched_samples
)

## 3. Oncoplots -------------------------------------------------------------

pdf(file.path(output_dir, "Oncoplot_SouthAsian.pdf"), width = 10, height = 8)
oncoplot(maf_filtered, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

pdf(file.path(output_dir, "Oncoplot_European.pdf"), width = 10, height = 8)
oncoplot(tcga_hnsc_eur_oscc, genes = comparison_genes, removeNonMutated = FALSE)
dev.off()

## 4. Gene frequency barplot ------------------------------------------------

gene_summary_sa <- getGeneSummary(maf_filtered) %>%
  mutate(Cohort = "South Asian")
gene_summary_eur <- getGeneSummary(tcga_hnsc_eur_oscc) %>%
  mutate(Cohort = "European")

combined_summary <- bind_rows(
  gene_summary_sa %>% filter(Hugo_Symbol %in% comparison_genes),
  gene_summary_eur %>% filter(Hugo_Symbol %in% comparison_genes)
)

p_freq <- ggplot(
  combined_summary,
  aes(x = reorder(Hugo_Symbol, -MutatedSamples),
      y = MutatedSamples,
      fill = Cohort)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Mutation Frequencies: South Asian vs European OSCC",
    x = "",
    y = "Mutated Samples"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_freq, "Mutation_Frequency_Comparison.pdf", width = 10, height = 6)

## 5. VAF distribution ------------------------------------------------------

df_vaf <- data.frame(
  VAF = c(maf_vaf$VAF, tcga_hnsc_eur_oscc@data$VAF),
  Cohort = c(
    rep("South Asian", length(maf_vaf$VAF)),
    rep("European", length(tcga_hnsc_eur_oscc@data$VAF))
  )
) %>%
  filter(!is.na(VAF))

p_vaf <- ggplot(df_vaf, aes(x = Cohort, y = VAF, fill = Cohort)) +
  geom_violin(alpha = 0.6, trim = FALSE, color = "black") +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.6, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  labs(
    title = "Variant Allele Frequency (VAF) Distribution",
    x = "",
    y = "VAF"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

save_plot(p_vaf, "VAF_Comparison.pdf", width = 7, height = 5)

## 6. Tumor mutational burden (TMB) ----------------------------------------

tmb_sa  <- tmb(maf_filtered)
tmb_eur <- tmb(tcga_hnsc_eur_oscc)

tmb_sa$Cohort  <- "South Asian"
tmb_eur$Cohort <- "European"

tmb_combined <- bind_rows(tmb_sa, tmb_eur)

p_tmb <- ggplot(tmb_combined, aes(x = Cohort, y = total_perMB, fill = Cohort)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "red", outlier.size = 2) +
  scale_fill_manual(values = c("South Asian" = "#66c2a5", "European" = "#fc8d62")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Tumor Mutational Burden (TMB) Comparison",
    x = "",
    y = "Mutations per MB"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_plot(p_tmb, "TMB_Comparison.pdf", width = 6, height = 5)

tmb_wilcox <- wilcox.test(total_perMB ~ Cohort, data = tmb_combined)
print(tmb_wilcox)

## 7. Fisher exact test + volcano plot -------------------------------------

get_fisher <- function(maf1, maf2) {
  g1 <- getGeneSummary(maf1)
  g2 <- getGeneSummary(maf2)
  
  total1 <- as.numeric(maf1@summary[maf1@summary$ID == "Samples", ]$summary)
  total2 <- as.numeric(maf2@summary[maf2@summary$ID == "Samples", ]$summary)
  
  all_genes <- union(g1$Hugo_Symbol, g2$Hugo_Symbol)
  fisher_table <- data.frame()
  
  for (g in all_genes) {
    a <- g1[g1$Hugo_Symbol == g, "MutatedSamples"]
    a <- if (length(a) == 0) 0 else as.numeric(a)
    
    c <- g2[g2$Hugo_Symbol == g, "MutatedSamples"]
    c <- if (length(c) == 0) 0 else as.numeric(c)
    
    b <- total1 - a
    d <- total2 - c
    
    vals <- c(a, b, c, d)
    if (any(is.na(vals)) || any(vals < 0) || any(!is.finite(vals))) next
    
    mat <- matrix(vals, nrow = 2, byrow = TRUE)
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
  fisher_table
}

fisher_table <- get_fisher(maf_filtered, tcga_hnsc_eur_oscc)
fisher_table <- fisher_table[is.finite(fisher_table$OR) & fisher_table$OR > 0, ]

fisher_table$logOR     <- log2(fisher_table$OR)
fisher_table$neglog10p <- -log10(fisher_table$p_adj)

write.csv(
  fisher_table,
  file.path(output_dir, "Fisher_Test_Results.csv"),
  row.names = FALSE
)

p_volcano <- ggplot(fisher_table, aes(x = logOR, y = neglog10p, color = p_adj < 0.05)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot – Differential Mutation Frequency",
    x = "log2(Odds Ratio)",
    y = "-log10(adjusted p-value)",
    color = "Significant"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_plot(p_volcano, "Volcano_Fisher.pdf", width = 7, height = 5)




# ============================================================
# 1. Define pathway gene sets
# ============================================================

wnt_genes  <- c("APC","CTNNB1","AXIN1","AXIN2","LRP5","LRP6","TCF7L2")
pi3k_genes <- c("PIK3CA","PIK3R1","PTEN","AKT1","AKT2","AKT3","MTOR","TSC1","TSC2")
ras_genes  <- c("KRAS","NRAS","HRAS","BRAF","MAPK1","MAP2K1")

# ============================================================
# 2. Robust Fisher function (version‑independent)
# ============================================================

get_fisher_clean <- function(maf1, maf2) {
  
  samples1 <- as.character(unique(maf1@data$Tumor_Sample_Barcode))
  samples2 <- as.character(unique(maf2@data$Tumor_Sample_Barcode))
  
  total1 <- length(samples1)
  total2 <- length(samples2)
  
  genes1 <- as.character(unique(maf1@data$Hugo_Symbol))
  genes2 <- as.character(unique(maf2@data$Hugo_Symbol))
  all_genes <- union(genes1, genes2)
  
  out <- data.frame()
  
  for (g in all_genes) {
    
    a <- length(unique(as.character(
      maf1@data$Tumor_Sample_Barcode[maf1@data$Hugo_Symbol == g]
    )))
    
    c <- length(unique(as.character(
      maf2@data$Tumor_Sample_Barcode[maf2@data$Hugo_Symbol == g]
    )))
    
    b <- total1 - a
    d <- total2 - c
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    if (any(mat == 0)) mat <- mat + 0.5
    
    ft <- fisher.test(mat)
    
    out <- rbind(
      out,
      data.frame(
        Gene = g,
        OR = as.numeric(ft$estimate),
        CI_low = ft$conf.int[1],
        CI_high = ft$conf.int[2],
        p = ft$p.value
      )
    )
  }
  
  out$p_adj <- p.adjust(out$p, method = "BH")
  out
}

# ============================================================
# 3. Run Fisher and annotate pathways + significance stars
# ============================================================

fisher_clean <- get_fisher_clean(maf_filtered, tcga_hnsc_eur_oscc)

forest_df <- fisher_clean %>%
  filter(!is.infinite(OR), OR > 0) %>%      # remove unusable ORs
  arrange(p_adj) %>%
  slice(1:25) %>%                           # top 25 most significant
  mutate(
    Pathway = case_when(
      Gene %in% wnt_genes  ~ "WNT",
      Gene %in% pi3k_genes ~ "PI3K",
      Gene %in% ras_genes  ~ "RAS",
      TRUE ~ "Other"
    ),
    Stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    GeneLabel = paste0(Gene, " ", Stars),
    GeneLabel = factor(GeneLabel, levels = rev(GeneLabel))
  )

# Pathway colors
pathway_colors <- c(
  "WNT"  = "#1b9e77",
  "PI3K" = "#d95f02",
  "RAS"  = "#7570b3",
  "Other" = "pink"
)

# ============================================================
# 4. Forest plot (clean, readable, colored by pathway)
# ============================================================

library(ggplot2)

p_forest <- ggplot(forest_df, aes(x = OR, y = GeneLabel, color = Pathway)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_log10() +
  scale_color_manual(values = pathway_colors) +
  theme_bw(base_size = 14) +
  labs(
    title = "European (n = 247) vs South Asian (n = 37)",
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 11),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "Forest_Plot_Euro_vs_SA_Pathways.pdf"),
  p_forest, width = 8, height = 6
)

save_plot(p_volcano, "Forest_Plot_Euro_vs_SA_Pathways.pdf", width = 8, height = 6)

# ============================================================
# Clean mutation frequency tables for scatter plot
# ============================================================

get_freq_table <- function(maf_obj) {
  df <- maf_obj@data
  
  # Count mutated samples per gene
  freq <- df %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(
      MutatedSamples = dplyr::n_distinct(Tumor_Sample_Barcode),
      .groups = "drop"
    )
  
  # Total samples in cohort
  total_samples <- length(unique(df$Tumor_Sample_Barcode))
  
  freq$Frequency <- freq$MutatedSamples / total_samples
  freq
}

sa_freq  <- get_freq_table(maf_filtered)
eur_freq <- get_freq_table(tcga_hnsc_eur_oscc)

# Merge cleanly
scatter_df <- dplyr::full_join(
  sa_freq %>% dplyr::rename(SA_Mut = MutatedSamples, SA_Freq = Frequency),
  eur_freq %>% dplyr::rename(EUR_Mut = MutatedSamples, EUR_Freq = Frequency),
  by = "Hugo_Symbol"
)

# Replace NAs with 0
scatter_df[is.na(scatter_df)] <- 0


# ============================================================
# Scatter Plot: Mutation Frequency (SA vs EUR)
# ============================================================

library(ggplot2)

p_scatter <- ggplot(scatter_df, aes(x = EUR_Freq, y = SA_Freq)) +
  geom_point(color = "#1f78b4", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  geom_text(
    aes(label = Hugo_Symbol),
    hjust = -0.1, vjust = 0.5, size = 3.2, check_overlap = TRUE
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Mutation Frequency: European vs South Asian",
    x = "European Mutation Frequency",
    y = "South Asian Mutation Frequency"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(output_dir, "Scatter_EUR_vs_SA.pdf"),
  p_scatter, width = 7.5, height = 6
)




## 8. Binary mutation matrix + PCA + UpSet ---------------------------------

# Build binary matrix for comparison_genes across both cohorts
build_binary_matrix <- function(maf_list, cohort_labels, genes) {
  genes <- toupper(trimws(genes))
  all_samples <- unlist(lapply(maf_list, function(m) unique(m@data$Tumor_Sample_Barcode)))
  all_samples <- unique(all_samples)
  
  bin_mat <- matrix(0, nrow = length(genes), ncol = length(all_samples),
                    dimnames = list(genes, all_samples))
  
  for (i in seq_along(maf_list)) {
    m <- maf_list[[i]]
    m@data$Hugo_Symbol <- toupper(trimws(m@data$Hugo_Symbol))
    for (g in genes) {
      mutated_samples <- unique(m@data$Tumor_Sample_Barcode[m@data$Hugo_Symbol == g])
      bin_mat[g, colnames(bin_mat) %in% mutated_samples] <- 1
    }
  }
  bin_mat
}

binary_matrix <- build_binary_matrix(
  maf_list = list(maf_filtered, tcga_hnsc_eur_oscc),
  cohort_labels = c("South Asian", "European"),
  genes = comparison_genes
)

# PCA
if (ncol(binary_matrix) >= 2) {
  pca_res <- prcomp(t(binary_matrix), scale. = TRUE)
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  
  sample_ids <- rownames(pca_res$x)
  sa_samples  <- unique(maf_filtered@data$Tumor_Sample_Barcode)
  eur_samples <- unique(tcga_hnsc_eur_oscc@data$Tumor_Sample_Barcode)
  
  pca_df <- data.frame(
    Sample = sample_ids,
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Cohort = case_when(
      Sample %in% sa_samples  ~ "South Asian",
      Sample %in% eur_samples ~ "European",
      TRUE ~ "Other"
    )
  )
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cohort, shape = Cohort)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 14) +
    labs(
      title = "PCA of Binary Mutation Matrix",
      x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "top"
    ) +
    geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 5)
  
  save_plot(p_pca, "PCA_BinaryMutation.pdf", width = 8, height = 6)
}


## Alternative version of PCA

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



# UpSet plot from binary matrix
upset_df <- as.data.frame(t(binary_matrix))
upset_df[upset_df > 1] <- 1

pdf(file.path(output_dir, "UpSet_BinaryMutation.pdf"), width = 10, height = 6)
upset(upset_df, nsets = min(10, ncol(upset_df)), nintersects = 30,
      order.by = "freq", mb.ratio = c(0.6, 0.4))
dev.off()

## 9. Variant classification pie charts s

plot_variant_pie <- function(maf_obj, cohort_name) {
  var_summary <- maf_obj@data %>%
    count(Variant_Classification) %>%
    mutate(Percent = n / sum(n) * 100)
  
  ggplot(var_summary, aes(x = "", y = Percent, fill = Variant_Classification)) +
    geom_col(color = "white") +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = paste0(round(Percent, 1), "%")),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    labs(
      title = paste("Variant Classification –", cohort_name),
      fill = "Mutation Type"
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

p_pie_sa  <- plot_variant_pie(maf_filtered, "South Asian")
p_pie_eur <- plot_variant_pie(tcga_hnsc_eur_oscc, "European")

save_plot(p_pie_sa,  "VariantClassification_SouthAsian.pdf", width = 7, height = 5)
save_plot(p_pie_eur, "VariantClassification_European.pdf", width = 7, height = 5)

## 10. Lollipop plots for key genes ----------------------------------------

genes_for_lollipop <- unique(c(wnt_genes, pi3k_genes, ras_genes, driver_genes))

for (gene in genes_for_lollipop) {
  if (gene %in% maf_filtered@data$Hugo_Symbol) {
    pdf(file.path(output_dir, paste0("Lollipop_", gene, "_SouthAsian.pdf")),
        width = 10, height = 4)
    lollipopPlot(
      maf_filtered,
      gene = gene,
      AACol = "HGVSp_Short",
      showMutationRate = TRUE
    )
    dev.off()
  }
}

## 11. Pathway mutation summary table --------------------------------------

pathway_summary <- data.frame(
  Pathway = c("WNT", "PI3K", "RAS"),
  Mutated_Samples_SouthAsian = sapply(
    list(wnt_genes, pi3k_genes, ras_genes),
    function(x) length(unique(
      maf_filtered@data$Tumor_Sample_Barcode[
        maf_filtered@data$Hugo_Symbol %in% x
      ]
    ))
  ),
  Mutated_Samples_European = sapply(
    list(wnt_genes, pi3k_genes, ras_genes),
    function(x) length(unique(
      tcga_hnsc_eur_oscc@data$Tumor_Sample_Barcode[
        tcga_hnsc_eur_oscc@data$Hugo_Symbol %in% x
      ]
    ))
  )
)

write.csv(
  pathway_summary,
  file.path(output_dir, "PathwayMutationTable.csv"),
  row.names = FALSE
)

## 12. KEGG enrichment (per cohort + comparison) ---------------------------

top_sa <- getGeneSummary(maf_filtered) %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 200) %>%
  pull(Hugo_Symbol)

top_tcga <- getGeneSummary(tcga_hnsc_eur_oscc) %>%
  arrange(desc(MutatedSamples)) %>%
  slice_head(n = 200) %>%
  pull(Hugo_Symbol)

gene2entrez <- bitr(
  unique(c(top_sa, top_tcga)),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_sa <- gene2entrez %>%
  filter(SYMBOL %in% top_sa) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

entrez_tcga <- gene2entrez %>%
  filter(SYMBOL %in% top_tcga) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

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

write.csv(as.data.frame(kegg_sa),
          file.path(output_dir, "KEGG_SA.csv"), row.names = FALSE)
write.csv(as.data.frame(kegg_tcga),
          file.path(output_dir, "KEGG_TCGA.csv"), row.names = FALSE)

# Clean SA KEGG plot (cancer pathways only)
cancer_keywords <- c("PI3K", "AKT", "WNT", "RAS", "MAPK", "p53", "Cell cycle", "Apoptosis")

kegg_sa_df <- as.data.frame(kegg_sa) %>%
  filter(grepl(paste(cancer_keywords, collapse = "|"),
               Description, ignore.case = TRUE)) %>%
  arrange(p.adjust) %>%
  mutate(Description = factor(Description, levels = Description))

if (nrow(kegg_sa_df) > 0) {
  p_kegg_sa <- ggplot(kegg_sa_df, aes(x = -log10(p.adjust), y = Description)) +
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
  
  save_plot(p_kegg_sa, "KEGG_Clean_SA_OSCC.pdf", width = 8, height = 6)
}

# Side-by-side KEGG comparison
kegg_sa_df2 <- as.data.frame(kegg_sa) %>%
  mutate(Cohort = "South Asian") %>%
  filter(grepl(paste(cancer_keywords, collapse = "|"),
               Description, ignore.case = TRUE))

kegg_tcga_df2 <- as.data.frame(kegg_tcga) %>%
  mutate(Cohort = "European") %>%
  filter(grepl(paste(cancer_keywords, collapse = "|"),
               Description, ignore.case = TRUE))

kegg_combined <- bind_rows(kegg_sa_df2, kegg_tcga_df2)

if (nrow(kegg_combined) > 0) {
  kegg_combined <- kegg_combined %>%
    arrange(Cohort, p.adjust) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p_kegg_comparison <- ggplot(
    kegg_combined,
    aes(x = -log10(p.adjust), y = Description, color = Cohort, size = Count)
  ) +
    geom_point(position = position_dodge(width = 0.6)) +
    scale_color_manual(values = c("South Asian" = "steelblue", "European" = "firebrick")) +
    theme_bw(base_size = 14) +
    labs(
      title = "KEGG Pathway Enrichment: South Asian vs European OSCC",
      x = "-log10(adjusted p-value)",
      y = "",
      color = "Cohort",
      size = "Mutated Genes"
    ) +
    theme(
      axis.text.y = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  save_plot(p_kegg_comparison, "KEGG_SideBySide_OSCC.pdf", width = 10, height = 6)
}


# ============================================
# Scatter plot: European vs South Asian mutation frequencies
# ============================================

# Robust frequency calculation
calc_freq <- function(maf, genes){
  total <- length(unique(maf@data$Tumor_Sample_Barcode))
  maf@data$Hugo_Symbol <- toupper(trimws(maf@data$Hugo_Symbol))
  genes <- toupper(trimws(genes))
  
  freq_df <- maf@data %>%
    dplyr::filter(Hugo_Symbol %in% genes) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(Samples = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
    dplyr::mutate(Frequency = Samples / total)
  
  # Add missing genes with frequency 0
  missing_genes <- setdiff(genes, freq_df$Hugo_Symbol)
  if (length(missing_genes) > 0) {
    freq_df <- dplyr::bind_rows(
      freq_df,
      data.frame(Hugo_Symbol = missing_genes,
                 Samples = 0,
                 Frequency = 0,
                 stringsAsFactors = FALSE)
    )
  }
  
  freq_df[order(-freq_df$Frequency), ]
}

# Calculate frequencies
south_freq <- calc_freq(maf_filtered, comparison_genes)
euro_freq  <- calc_freq(tcga_hnsc_eur_oscc, comparison_genes)

# Base R merge to avoid dplyr name repair
m <- merge(
  south_freq[, c("Hugo_Symbol", "Frequency")],
  euro_freq[,  c("Hugo_Symbol", "Frequency")],
  by = "Hugo_Symbol",
  all = TRUE,
  suffixes = c("_South", "_Europe")
)

# Build clean freq_df with unique column names
freq_df <- data.frame(
  Hugo_Symbol = m$Hugo_Symbol,
  South       = ifelse(is.na(m$Frequency_South),  0, m$Frequency_South),
  Europe      = ifelse(is.na(m$Frequency_Europe), 0, m$Frequency_Europe),
  stringsAsFactors = FALSE
)

# Assign pathway categories
freq_df$Pathway <- "Other"
freq_df$Pathway[freq_df$Hugo_Symbol %in% wnt_genes]  <- "WNT"
freq_df$Pathway[freq_df$Hugo_Symbol %in% pi3k_genes] <- "PI3K"
freq_df$Pathway[freq_df$Hugo_Symbol %in% ras_genes]  <- "RAS"

# Scatter plot
library(ggrepel)

p_scatter <- ggplot(freq_df, aes(x = Europe, y = South, color = Pathway, label = Hugo_Symbol)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c(
    "WNT"   = "#1f78b4",
    "PI3K"  = "#33a02c",
    "RAS"   = "#e31a1c",
    "Other" = "pink"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Mutation Frequency Comparison: South Asian vs European OSCC",
    x = "European OSCC Frequency",
    y = "South Asian OSCC Frequency",
    color = "Pathway"
  ) +
  xlim(0, 1) + ylim(0, 1)

ggsave(file.path(output_dir, "Comparison_Frequency_Scatter_Pretty.pdf"),
       p_scatter, width = 9, height = 6)

write.csv(freq_df, file.path(output_dir, "Mutation_Frequency_Table.csv"), row.names = FALSE)



# ============================================================
# 1. Compute mutation-type proportions for each cohort
# ============================================================

get_mutation_type_table <- function(maf_obj, cohort_name) {
  df <- maf_obj@data
  
  type_counts <- df %>%
    dplyr::group_by(Variant_Classification) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      prop = n / sum(n),
      Cohort = cohort_name
    )
  
  type_counts
}

sa_types  <- get_mutation_type_table(maf_filtered, "South Asian")
eur_types <- get_mutation_type_table(tcga_hnsc_eur_oscc, "European")

combined_types <- dplyr::bind_rows(sa_types, eur_types)

mut_colors <- c(
  "Missense_Mutation"      = "#1b9e77",
  "Nonsense_Mutation"      = "#d95f02",
  "Frame_Shift_Del"        = "#7570b3",
  "Frame_Shift_Ins"        = "#e7298a",
  "In_Frame_Del"           = "#66a61e",
  "In_Frame_Ins"           = "#e6ab02",
  "Splice_Site"            = "#a6761d",
  "Translation_Start_Site" = "#666666",
  "Nonstop_Mutation"       = "#1f78b4",
  "Silent"                 = "#b2df8a"
)

library(ggplot2)

p_stack <- ggplot(
  combined_types,
  aes(x = Cohort, y = prop, fill = Variant_Classification)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = mut_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 14) +
  labs(
    title = "Mutation Type Distribution in OSCC: South Asian vs European",
    x = "",
    y = "Proportion of Mutations"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "StackedBar_OSCC_SA_vs_EUR.pdf"),
  p_stack, width = 8, height = 6
)




## 14. Save core objects & session info ------------------------------------

save(
  maf_filtered,
  tcga_hnsc_eur_oscc,
  tmb_combined,
  fisher_table,
  kegg_sa,
  kegg_tcga,
  binary_matrix,
  file = file.path(output_dir, "all_objects.RData")
)



# ============================================
# VAF Comparison – South Asian vs European
# ============================================

# Extract VAFs safely
vaf_sa   <- maf_vaf$VAF
vaf_eur  <- tcga_hnsc_eur_oscc@data$VAF

# Build clean dataframe
df_vaf <- data.frame(
  VAF = c(vaf_sa, vaf_eur),
  Cohort = c(
    rep("South Asian", length(vaf_sa)),
    rep("European",   length(vaf_eur))
  ),
  stringsAsFactors = FALSE
)

# Remove NA values
df_vaf <- df_vaf[!is.na(df_vaf$VAF), ]

# Compute VAF for TCGA European OSCC if missing
if (!"VAF" %in% colnames(tcga_hnsc_eur_oscc@data)) {
  if (all(c("t_alt_count", "t_ref_count") %in% colnames(tcga_hnsc_eur_oscc@data))) {
    
    tcga_hnsc_eur_oscc@data$VAF <- with(
      tcga_hnsc_eur_oscc@data,
      ifelse(
        (t_alt_count + t_ref_count) > 0,
        t_alt_count / (t_alt_count + t_ref_count),
        NA
      )
    )
    
  } else {
    stop("European TCGA MAF does not contain t_alt_count and t_ref_count — cannot compute VAF.")
  }
}


# Pretty VAF plot
library(ggplot2)

p_vaf <- ggplot(df_vaf, aes(x = Cohort, y = VAF, fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.5, color = "black") +
  stat_summary(fun = median, geom = "point",
               shape = 23, size = 3.5, fill = "red", color = "black") +
  scale_fill_manual(values = c(
    "South Asian" = "#66c2a5",
    "European"    = "#fc8d62"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 15) +
  labs(
    title = "Variant Allele Frequency (VAF) Distribution",
    x = "",
    y = "Variant Allele Frequency (VAF)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text = element_text(size = 13),
    legend.position = "none"
  )

# Save
ggsave(
  file.path(output_dir, "VAF_Comparison_Pretty.pdf"),
  p_vaf, width = 7.5, height = 5.5
)

cat("\n===== PIPELINE COMPLETED SUCCESSFULLY =====\n\n")
sessionInfo()

library(ggplot2)

p_stack <- ggplot(
  combined_types,
  aes(x = Cohort, y = prop, fill = Variant_Classification)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = mut_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 14) +
  labs(
    title = "Mutation Type Distribution in OSCC: South Asian vs European",
    x = "",
    y = "Proportion of Mutations"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "StackedBar_OSCC_SA_vs_EUR.pdf"),
  p_stack, width = 8, height = 6
)






# ============================================================
# Oncoplot for top mutated genes (SA vs EUR)
# ============================================================

library(maftools)

# Top 15 genes across both cohorts
top_genes <- fisher_clean %>%
  arrange(p_adj) %>%
  slice(1:15) %>%
  pull(Gene)

# South Asian oncoplot
oncoplot(
  maf = maf_filtered,
  genes = top_genes,
  draw_titv = FALSE,
  sortByAnnotation = TRUE,
  annotationColor = NULL,
  removeNonMutated = TRUE,
  titleText = "South Asian OSCC – Top Mutated Genes"
)

# European oncoplot
oncoplot(
  maf = tcga_hnsc_eur_oscc,
  genes = top_genes,
  draw_titv = FALSE,
  sortByAnnotation = TRUE,
  annotationColor = NULL,
  removeNonMutated = TRUE,
  titleText = "European OSCC – Top Mutated Genes"
)


install.packages("NMF")
BiocManager::install("maftools")

install.packages("NMF", type = "source")

library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38)


# ============================================================
# 1. Build trinucleotide matrix
# ============================================================

maf_filtered@data$Chromosome <- paste0("chr", maf_filtered@data$Chromosome)

sa_tnm <- trinucleotideMatrix(
  maf = maf_filtered,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
)

tcga_hnsc_eur_oscc@data$Chromosome <- paste0("chr", tcga_hnsc_eur_oscc@data$Chromosome)
eur_tnm <- trinucleotideMatrix(
  maf = tcga_hnsc_eur_oscc,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
)

# ============================================================
# 2. Extract signatures (NMF)
# ============================================================

library(maftools)
library(NMF)

# sa_tnm is already created and valid

nmf_mat <- sa_tnm$nmf_matrix   # extract the matrix safely

sa_sig <- extractSignatures(
  mat = nmf_mat,
  n = 5
)

plotSignatures(sa_sig, title_size = 1.2)


# ============================================================
# 3. Compare to COSMIC SBS reference signatures
# ============================================================

sa_cosmic <- compareSignatures(
  sa_sig$signatures,
  sig_db = "SBS"
)

eur_cosmic <- compareSignatures(
  eur_sig$signatures,
  sig_db = "SBS"
)

# ============================================================
# 4. Plot similarity heatmaps
# ============================================================

plot(sa_cosmic)
plot(eur_cosmic)


library(maftools)

# Step 1: build trinucleotide matrix (you already did this)
sa_tnm <- trinucleotideMatrix(
  maf = maf_filtered,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
)

# Step 2: estimate signatures (this replaces extractSignatures)
sa_sig <- estimateSignatures(
  mat = sa_tnm$nmf_matrix,
  n = 5
)

# Step 3: plot signatures
plotSignatures(sa_sig)









# Remove Duplicates

fix_all_maf_colnames <- function(maf_obj) {
  # Fix main data slot
  cn <- colnames(maf_obj@data)
  colnames(maf_obj@data) <- make.unique(cn, sep = "_")
  
  # Fix clinical data slot (if exists)
  if (!is.null(maf_obj@clinical.data)) {
    cn2 <- colnames(maf_obj@clinical.data)
    colnames(maf_obj@clinical.data) <- make.unique(cn2, sep = "_")
  }
  
  # Fix summary slot (rarely needed but safe)
  if (!is.null(maf_obj@summary)) {
    cn3 <- colnames(maf_obj@summary)
    colnames(maf_obj@summary) <- make.unique(cn3, sep = "_")
  }
  
  maf_obj
}

maf_filtered <- fix_all_maf_colnames(maf_filtered)
tcga_hnsc_eur_oscc <- fix_all_maf_colnames(tcga_hnsc_eur_oscc)


# Statistical tests

vars_to_test <- c("TMB", "MedianVAF", "MutationCount", "Age")

# Convert to numeric safely
df[, vars_to_test] <- lapply(df[, vars_to_test], function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
})

# Remove rows with NA in any test variable
df_clean <- df[complete.cases(df[, vars_to_test]), ]


library(car)

compare_two_cohorts_safe <- function(df, value_col, group_col = "Cohort") {
  
  g1 <- unique(df[[group_col]])[1]
  g2 <- unique(df[[group_col]])[2]
  
  x <- df[df[[group_col]] == g1, value_col]
  y <- df[df[[group_col]] == g2, value_col]
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  n1 <- length(x)
  n2 <- length(y)
  
  # If either group has < 2 samples → cannot test
  if (n1 < 2 || n2 < 2) {
    return(list(
      method = "Insufficient sample size",
      p_value = NA,
      statistic = NA,
      n_group1 = n1,
      n_group2 = n2
    ))
  }
  
  # If either group has < 3 samples → skip Shapiro, use Wilcoxon
  if (n1 < 3 || n2 < 3) {
    test <- wilcox.test(x, y)
    return(list(
      method = "Wilcoxon (small sample fallback)",
      p_value = test$p.value,
      statistic = test$statistic,
      n_group1 = n1,
      n_group2 = n2
    ))
  }
  
  # Normality tests
  shapiro_x <- shapiro.test(x)
  shapiro_y <- shapiro.test(y)
  
  normal_x <- shapiro_x$p.value > 0.05
  normal_y <- shapiro_y$p.value > 0.05
  
  # Variance test only if both normal
  if (normal_x && normal_y) {
    lev <- car::leveneTest(df[[value_col]] ~ df[[group_col]])
    equal_var <- lev$`Pr(>F)`[1] > 0.05
  } else {
    equal_var <- FALSE
  }
  
  # Choose test
  if (normal_x && normal_y) {
    test <- t.test(x, y, var.equal = equal_var)
    method <- ifelse(equal_var, "Student t-test", "Welch t-test")
  } else {
    test <- wilcox.test(x, y)
    method <- "Wilcoxon rank-sum test"
  }
  
  return(list(
    method = method,
    p_value = test$p.value,
    statistic = test$statistic,
    shapiro_x = shapiro_x$p.value,
    shapiro_y = shapiro_y$p.value,
    equal_variance = equal_var,
    n_group1 = n1,
    n_group2 = n2
  ))
}


run_all_tests_safe <- function(df, vars, group_col = "Cohort") {
  results <- lapply(vars, function(v) {
    out <- compare_two_cohorts_safe(df, value_col = v, group_col = group_col)
    data.frame(
      Variable = v,
      Method = out$method,
      P_value = out$p_value,
      Statistic = out$statistic,
      N_group1 = out$n_group1,
      N_group2 = out$n_group2
    )
  })
  do.call(rbind, results)
}


vars_to_test <- c("VAF")
results <- run_all_tests_safe(df, vars_to_test)
results



# ============================================================
# 2. Robust Fisher function (version-independent)
# ============================================================

get_fisher_clean <- function(maf1, maf2) {
  
  samples1 <- as.character(unique(maf1@data$Tumor_Sample_Barcode))
  samples2 <- as.character(unique(maf2@data$Tumor_Sample_Barcode))
  
  total1 <- length(samples1)
  total2 <- length(samples2)
  
  genes1 <- as.character(unique(maf1@data$Hugo_Symbol))
  genes2 <- as.character(unique(maf2@data$Hugo_Symbol))
  all_genes <- union(genes1, genes2)
  
  out <- data.frame()
  
  for (g in all_genes) {
    
    a <- length(unique(as.character(
      maf1@data$Tumor_Sample_Barcode[maf1@data$Hugo_Symbol == g]
    )))
    
    c <- length(unique(as.character(
      maf2@data$Tumor_Sample_Barcode[maf2@data$Hugo_Symbol == g]
    )))
    
    b <- total1 - a
    d <- total2 - c
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    if (any(mat == 0)) mat <- mat + 0.5
    
    ft <- fisher.test(mat)
    
    out <- rbind(
      out,
      data.frame(
        Gene = g,
        OR = as.numeric(ft$estimate),
        CI_low = ft$conf.int[1],
        CI_high = ft$conf.int[2],
        p = ft$p.value
      )
    )
  }
  
  out$p_adj <- p.adjust(out$p, method = "BH")
  out
}

# ============================================================
# 3. Run Fisher and annotate pathways BEFORE slicing
# ============================================================

fisher_clean <- get_fisher_clean(maf_filtered, tcga_hnsc_eur_oscc)

# Force plain tibble to avoid hidden Bioconductor attributes
fisher_clean <- fisher_clean |> tibble::as_tibble()

# Annotate first, THEN slice
annotated_df <- fisher_clean %>%
  filter(!is.infinite(OR), OR > 0) %>%
  arrange(p_adj) %>%
  mutate(
    Pathway = case_when(
      Gene %in% wnt_genes  ~ "WNT",
      Gene %in% pi3k_genes ~ "PI3K",
      Gene %in% ras_genes  ~ "RAS",
      TRUE ~ "Other"
    ),
    Stars = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    GeneLabel = paste0(Gene, " ", Stars)
  )

annotated_df <- annotated_df %>%
  mutate(
    Gene   = as.character(Gene),
    OR     = as.numeric(OR),
    CI_low = as.numeric(CI_low),
    CI_high= as.numeric(CI_high),
    p      = as.numeric(p),
    p_adj  = as.numeric(p_adj),
    Pathway = as.character(Pathway),
    Stars   = as.character(Stars),
    GeneLabel = as.character(GeneLabel)
  )

# FULL hard reset of all internal attributes:
annotated_df <- data.frame(annotated_df, stringsAsFactors = FALSE)


# Slice AFTER all vectors have proper length
forest_df <- annotated_df %>%
  .[seq_len(min(25, nrow(.))), ]
  mutate(GeneLabel = factor(GeneLabel, levels = rev(GeneLabel)))

# ============================================================
# Pathway colors
# ============================================================

pathway_colors <- c(
  "WNT"  = "#1b9e77",
  "PI3K" = "#d95f02",
  "RAS"  = "#7570b3",
  "Other" = "pink"
)

# ============================================================
# 4. Forest plot
# ============================================================

library(ggplot2)

p_forest <- ggplot(forest_df, aes(x = OR, y = GeneLabel, color = Pathway)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_log10() +
  scale_color_manual(values = pathway_colors) +
  theme_bw(base_size = 14) +
  labs(
    title = "European (n = 436) vs South Asian (n = 37)",
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 11),
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "Forest_Plot_Euro_vs_SA_Pathways.pdf"),
  p_forest, width = 8, height = 6
)
