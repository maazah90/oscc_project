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
output_dir <- "maf_results"
dir.create(output_dir, showWarnings = FALSE)

# --------------------------------
# Mutation burden per sample
# --------------------------------
south_mb <- getSampleSummary(maf_filtered) %>%
  select(Tumor_Sample_Barcode, total) %>%
  mutate(Cohort = "South Asian")

euro_mb <- getSampleSummary(tcga_hnsc_eur) %>%
  select(Tumor_Sample_Barcode, total) %>%
  mutate(Cohort = "European")

mb_df <- bind_rows(south_mb, euro_mb)

# --------------------------------
# Plot (publication style)
# --------------------------------
p <- ggplot(mb_df, aes(x = Cohort, y = total, fill = Cohort)) +
  
  geom_violin(trim = FALSE, alpha = 0.5) +
  
  geom_boxplot(width = 0.2, outlier.color = "red", fill = "white") +
  
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Mutation Burden Comparison",
    y = "Mutations per Sample",
    x = ""
  )

ggsave(file.path(output_dir, "Comparison_MutationBurden.pdf"),
       p, width = 6, height = 5)


get_pathway_freq <- function(maf, pathway_genes, name){
  
  data.frame(
    Pathway = name,
    Mutated_Samples = length(unique(
      maf@data$Tumor_Sample_Barcode[
        maf@data$Hugo_Symbol %in% pathway_genes
      ]
    )),
    Total_Samples = length(unique(maf@data$Tumor_Sample_Barcode))
  ) %>%
    mutate(Frequency = Mutated_Samples / Total_Samples)
}

# South Asian
south_pathways <- bind_rows(
  get_pathway_freq(maf_filtered, wnt_genes, "WNT"),
  get_pathway_freq(maf_filtered, pi3k_genes, "PI3K"),
  get_pathway_freq(maf_filtered, ras_genes, "RAS")
) %>% mutate(Cohort = "South Asian")

# European
euro_pathways <- bind_rows(
  get_pathway_freq(tcga_hnsc_eur, wnt_genes, "WNT"),
  get_pathway_freq(tcga_hnsc_eur, pi3k_genes, "PI3K"),
  get_pathway_freq(tcga_hnsc_eur, ras_genes, "RAS")
) %>% mutate(Cohort = "European")

pathway_df <- bind_rows(south_pathways, euro_pathways)

# Plot
p <- ggplot(pathway_df, aes(x = Pathway, y = Frequency, fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic(base_size = 14) +
  labs(
    title = "Pathway Alteration Frequency",
    y = "Proportion of Samples",
    x = ""
  )

ggsave(file.path(output_dir, "Comparison_Pathway_Frequency.pdf"),
       p, width = 6, height = 5)


# --------------------------------
# Compute frequencies
# --------------------------------
# --------------------------------
# FIX TCGA duplicate columns
# --------------------------------

library(data.table)

# --------------------------------
# FIX TCGA duplicate columns (SAFE)
# --------------------------------
tcga_data <- tcga_hnsc_eur@data

# Remove duplicate columns
tcga_data <- tcga_data[, !duplicated(colnames(tcga_data)), with = FALSE]

# Ensure it's still a data.table
setDT(tcga_data)

# Put back into MAF
tcga_hnsc_eur@data <- tcga_data

calc_freq <- function(maf, genes){
  
  total <- length(unique(maf@data$Tumor_Sample_Barcode))
  
  maf@data %>%
    filter(Hugo_Symbol %in% genes) %>%
    group_by(Hugo_Symbol) %>%
    summarise(Samples = n_distinct(Tumor_Sample_Barcode)) %>%
    mutate(Frequency = Samples / total)
}

south_freq <- calc_freq(maf_filtered, comparison_genes) %>%
  rename(South = Frequency)

euro_freq <- calc_freq(tcga_hnsc_eur, comparison_genes) %>%
  rename(Europe = Frequency)

freq_df <- full_join(south_freq, euro_freq, by = "Hugo_Symbol") %>%
  replace_na(list(South = 0, Europe = 0))

# --------------------------------
# Scatter plot
# --------------------------------
p <- ggplot(freq_df, aes(x = Europe, y = South)) +
  
  geom_point(size = 3, alpha = 0.7) +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  
  geom_text(aes(label = Hugo_Symbol), size = 3, vjust = -0.5) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Gene Mutation Frequency Comparison",
    x = "European Frequency",
    y = "South Asian Frequency"
  )

ggsave(file.path(output_dir, "Comparison_GeneScatter.pdf"),
       p, width = 6, height = 6)


df <- data.frame(
  VAF = c(vaf1, vaf2),
  Cohort = c(
    rep("South Asian", length(vaf1)),
    rep("European", length(vaf2))
  )
)

p <- ggplot(df, aes(x = Cohort, y = VAF, fill = Cohort)) +
  
  geom_violin(trim = FALSE, alpha = 0.5) +
  
  geom_boxplot(width = 0.2, fill = "white") +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Variant Allele Frequency Comparison",
    y = "VAF",
    x = ""
  )

ggsave(file.path(output_dir, "Comparison_VAF_Boxplot.pdf"),
       p, width = 6, height = 5)





install.packages("ggrepel")

library(ggrepel)

# --------------------------------
# Add pathway annotation
# --------------------------------
freq_df <- freq_df %>%
  mutate(
    Pathway = case_when(
      Hugo_Symbol %in% wnt_genes ~ "WNT",
      Hugo_Symbol %in% pi3k_genes ~ "PI3K",
      Hugo_Symbol %in% ras_genes ~ "RAS",
      TRUE ~ "Other"
    )
  )

# --------------------------------
# Select genes to label (important ones only)
# --------------------------------
genes_to_label <- freq_df %>%
  filter(
    Hugo_Symbol %in% c("TP53","PIK3CA","CDKN2A","NOTCH1","FAT1","CTNNB1") |
      abs(South - Europe) > 0.1   # big differences
  )

# --------------------------------
# Plot
# --------------------------------
p <- ggplot(freq_df, aes(x = Europe, y = South, color = Pathway)) +
  
  geom_point(size = 3, alpha = 0.8) +
  
  # diagonal reference
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  
  # 🔥 smart labels (no overlap)
  geom_text_repel(
    data = genes_to_label,
    aes(label = Hugo_Symbol),
    size = 3,
    max.overlaps = 50
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Gene Mutation Frequency Comparison",
    x = "European Frequency",
    y = "South Asian Frequency",
    color = "Pathway"
  )

ggsave(file.path(output_dir, "Comparison_GeneScatter_Improved.pdf"),
       p, width = 7, height = 6)



variant_counts <- maf_filtered@data %>%
  count(Variant_Classification) %>%
  arrange(desc(n))

p_variant <- ggplot(variant_counts,
                    aes(x = reorder(Variant_Classification, n),
                        y = n,
                        fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Variant Classification Distribution",
    x = "",
    y = "Number of Mutations"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(output_dir, "Fig_Variant_Barplot.pdf"),
       p_variant, width = 7, height = 5)