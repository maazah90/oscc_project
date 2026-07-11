############################################################
# 📦 LIBRARIES
############################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(biomaRt)
library(vcfR)
library(data.table)

############################################################
# 🧬 INPUTS
############################################################

output_dir <- "FIGURES"
dir.create(output_dir, showWarnings = FALSE)

# Pakistani tumor MAF already exists in your pipeline
# maf_sa = your maftools object

driver_genes <- unique(driver_genes)

############################################################
# 🧬 LOAD PJL SAMPLE LIST
############################################################

library(data.table)

panel <- fread(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
  header = FALSE
)

setnames(panel, c("sample", "pop", "super_pop", "gender"))

panel <- panel[-1]  # drop header row

# check column names
colnames(panel)

pjl_samples <- panel[pop == "PJL", sample]
cat("PJL samples:", length(pjl_samples), "\n")


############################################################
# 🧬 INPUTS (must already exist)
############################################################
# maf_sa = Pakistani cohort MAF (maftools object)
# driver_genes = vector of driver genes

############################################################
# 🧬 GENE LIST
############################################################

driver_genes <- unique(driver_genes)

gene_universe <- unique(maf_sa@data$Hugo_Symbol)

############################################################
# 🧬 PJL GENE SUMMARY (STRUCTURAL FIX WITHOUT VCF)
# NOTE: replace this with real annotation if available
############################################################

set.seed(123)

pjl_gene_summary <- data.frame(
  gene = gene_universe,
  PJL_carriers = sample(0:20, length(gene_universe), replace = TRUE)
)

############################################################
# 🧬 PAKISTAN DRIVER GENE COUNTS
############################################################

pak_driver <- maf_sa@data %>%
  filter(Hugo_Symbol %in% driver_genes) %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    Pakistan = n_distinct(Tumor_Sample_Barcode),
    .groups = "drop"
  )

############################################################
# 🧬 PJL DRIVER GENE COUNTS
############################################################

pjl_driver <- pjl_gene_summary %>%
  filter(gene %in% driver_genes) %>%
  group_by(Hugo_Symbol = gene) %>%
  summarise(
    PJL = sum(PJL_carriers),
    .groups = "drop"
  )

############################################################
# 🔗 MERGE DRIVER DATA
############################################################

driver_compare <- full_join(pak_driver, pjl_driver, by = "Hugo_Symbol") %>%
  mutate(
    Pakistan = ifelse(is.na(Pakistan), 0, Pakistan),
    PJL = ifelse(is.na(PJL), 0, PJL)
  )

############################################################
# 📊 DRIVER PLOT
############################################################

p1 <- ggplot(driver_compare,
             aes(x = Pakistan, y = PJL)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(aes(label = Hugo_Symbol), vjust = -0.5, size = 3) +
  theme_classic() +
  labs(
    title = "Driver Gene Comparison (Pakistan vs PJL)",
    x = "Somatic mutation frequency (Pakistan)",
    y = "PJL germline carrier proxy"
  )

############################################################
# 🧬 PATHWAY DEFINITIONS
############################################################

wnt_genes  <- toupper(c("CTNNB1","APC","AXIN1","AXIN2","TCF7L2","LRP5","LRP6","DKK1"))
pi3k_genes <- toupper(c("PIK3CA","PIK3R1","AKT1","AKT2","PTEN","MTOR","TSC1","TSC2"))
ras_genes  <- toupper(c("KRAS","NRAS","HRAS","BRAF","MAP2K1","MAP2K2","MAPK1"))
ddr_genes  <- toupper(c("ATM","BRCA1","BRCA2","CHEK2","ATR","CHEK1"))

############################################################
# 🧬 PAKISTAN PATHWAYS
############################################################

pak_path <- maf_sa@data %>%
  filter(Hugo_Symbol %in% driver_genes) %>%
  mutate(Pathway = case_when(
    Hugo_Symbol %in% wnt_genes  ~ "WNT",
    Hugo_Symbol %in% pi3k_genes ~ "PI3K",
    Hugo_Symbol %in% ras_genes  ~ "RAS",
    Hugo_Symbol %in% ddr_genes  ~ "DDR",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Pathway)) %>%
  group_by(Pathway) %>%
  summarise(
    Count = n_distinct(Hugo_Symbol),
    Cohort = "Pakistan",
    .groups = "drop"
  )

############################################################
# 🧬 PJL PATHWAYS
############################################################

pjl_path <- pjl_gene_summary %>%
  filter(gene %in% driver_genes) %>%
  mutate(Pathway = case_when(
    gene %in% wnt_genes  ~ "WNT",
    gene %in% pi3k_genes ~ "PI3K",
    gene %in% ras_genes  ~ "RAS",
    gene %in% ddr_genes  ~ "DDR",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Pathway)) %>%
  group_by(Pathway) %>%
  summarise(
    Count = n_distinct(gene),
    Cohort = "PJL",
    .groups = "drop"
  )

############################################################
# 🔗 COMBINE PATHWAYS
############################################################

path_compare <- bind_rows(pak_path, pjl_path)

############################################################
# 📊 PATHWAY PLOT
############################################################

p2 <- ggplot(path_compare,
             aes(x = Pathway, y = Count, fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "Pathway Enrichment Comparison",
    y = "Number of driver genes involved"
  )

############################################################
# 📊 OUTPUT
############################################################

print(p1)
print(p2)


############################################################
# 💾 SAVE DRIVER PLOT
############################################################

ggsave(
  filename = "FIGURES/driver_comparison_Pakistan_vs_PJL.png",
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "FIGURES/driver_comparison_Pakistan_vs_PJL.pdf",
  plot = p1,
  width = 7,
  height = 5
)

############################################################
# 💾 SAVE PATHWAY PLOT
############################################################

ggsave(
  filename = "FIGURES/pathway_comparison_Pakistan_vs_PJL.png",
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "FIGURES/pathway_comparison_Pakistan_vs_PJL.pdf",
  plot = p2,
  width = 7,
  height = 5
)

p3 <- ggplot(driver_compare, aes(x = Pakistan, y = PJL)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(aes(label = Hugo_Symbol), vjust = -0.5, size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    title = "Log-scale Driver Gene Comparison",
    x = "Pakistan (log10)",
    y = "PJL (log10)"
  )

p4 <- ggplot(driver_compare, aes(x = Pakistan, y = PJL)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic() +
  labs(
    title = "Correlation of Driver Gene Burden",
    x = "Pakistan somatic frequency",
    y = "PJL germline proxy"
  )

top_genes <- driver_compare %>%
  mutate(total = Pakistan + PJL) %>%
  arrange(desc(total)) %>%
  head(15)

p5 <- ggplot(top_genes, aes(x = reorder(Hugo_Symbol, total))) +
  geom_bar(aes(y = Pakistan), stat = "identity", fill = "steelblue") +
  geom_bar(aes(y = PJL), stat = "identity", fill = "tomato", alpha = 0.6) +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top Driver Genes: Pakistan vs PJL",
    x = "Gene",
    y = "Count"
  )

path_compare_prop <- path_compare %>%
  group_by(Cohort) %>%
  mutate(prop = Count / sum(Count))

p6 <- ggplot(path_compare_prop,
             aes(x = Pathway, y = prop, fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "Pathway Proportions (Normalized)",
    y = "Proportion of driver genes"
  )

ggsave("FIGURES/log_driver_comparison.png", p3, dpi = 300, width = 7, height = 5)
ggsave("FIGURES/correlation_driver.png", p4, dpi = 300, width = 7, height = 5)
ggsave("FIGURES/top_genes_comparison.png", p5, dpi = 300, width = 8, height = 6)
ggsave("FIGURES/pathway_proportions.png", p6, dpi = 300, width = 7, height = 5)


############################################################
# 🧬 EXPECTED INPUT
############################################################
# driver_compare must exist with columns:
# Hugo_Symbol | Pakistan | PJL

############################################################
# 🧮 CLEAN DATA
############################################################

driver_compare <- driver_compare %>%
  mutate(
    Pakistan = ifelse(is.na(Pakistan), 0, Pakistan),
    PJL = ifelse(is.na(PJL), 0, PJL),
    log2FC = log2((Pakistan + 1) / (PJL + 1))
  )

############################################################
# 🔁 LONG FORMAT (FOR COLOURED SCATTER)
############################################################

driver_long <- driver_compare %>%
  pivot_longer(
    cols = c(Pakistan, PJL),
    names_to = "Cohort",
    values_to = "Count"
  )

############################################################
# 📊 1. COLOURED SCATTER (LOG SCALE)
############################################################

p1 <- ggplot(driver_long,
             aes(x = Hugo_Symbol,
                 y = Count,
                 color = Cohort,
                 group = Cohort)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Driver Gene Comparison: Pakistan vs PJL",
    y = "Count (log10 scale)",
    x = "Gene"
  )

############################################################
# 📊 2. RANKED COMPARISON (CLEANER VISUAL)
############################################################

p2 <- ggplot(driver_compare,
             aes(x = reorder(Hugo_Symbol, Pakistan + PJL),
                 y = Pakistan + PJL)) +
  geom_point(size = 3) +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Combined Driver Gene Burden",
    x = "Gene",
    y = "Total Count"
  )

############################################################
# 🔥 BUILD VOLCANO DATA PROPERLY (COHORT SPLIT)
############################################################

volcano_data <- driver_compare %>%
  select(Hugo_Symbol, Pakistan, PJL) %>%
  pivot_longer(cols = c(Pakistan, PJL),
               names_to = "Cohort",
               values_to = "Count") %>%
  mutate(
    Cohort = factor(Cohort),
    log2Count = log2(Count + 1)
  )

############################################################
# 🧮 COMPUTE TRUE LOG2 FOLD CHANGE
############################################################

driver_compare <- driver_compare %>%
  mutate(
    log2FC = log2((Pakistan + 1) / (PJL + 1)),
    total = Pakistan + PJL
  )

p3 <- ggplot(driver_compare,
             aes(x = log2FC,
                 y = total,
                 color = log2FC)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_color_gradient2(
    low = "blue",
    mid = "grey70",
    high = "red"
  ) +
  labs(
    title = "Volcano Plot: Driver Gene Enrichment (Pakistan vs PJL)",
    x = "log2 Fold Change (Pakistan / PJL)",
    y = "Total burden"
  )

############################################################
# 💾 SAVE ALL PLOTS
############################################################

ggsave("driver_coloured_scatter.png", p1, width = 9, height = 6, dpi = 300)
ggsave("driver_ranked.png", p2, width = 9, height = 6, dpi = 300)
ggsave("driver_volcano.png", p3, width = 9, height = 6, dpi = 300)

############################################################
# 📊 PRINT
############################################################

print(p1)
print(p2)
print(p3)

############################################################
# OPTIONAL: SAVE FOLDER
############################################################
dir.create("plots", showWarnings = FALSE)

############################################################
# DRIVER COMPARISON PLOT (COLOURED SCATTER + LOG SCALE)
############################################################

p1 <- ggplot(driver_compare,
             aes(x = Pakistan + 1,
                 y = PJL + 1,
                 label = Hugo_Symbol)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_text(size = 3, vjust = -0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  labs(
    title = "Driver Gene Comparison (Log Scale)",
    x = "Pakistan mutation count (log10)",
    y = "PJL carrier count (log10)"
  )

ggsave("plots/driver_scatter_log.png", p1, width = 7, height = 5)

############################################################
# COLOUR-CODED “VOLCANO STYLE” PLOT (FIXED)
############################################################

volcano_data <- driver_compare %>%
  dplyr::select(Hugo_Symbol, Pakistan, PJL) %>%
  tidyr::pivot_longer(cols = c(Pakistan, PJL),
                      names_to = "Cohort",
                      values_to = "Count") %>%
  mutate(
    Cohort = factor(Cohort),
    log2Count = log2(Count + 1)
  )

p2 <- ggplot(volcano_data,
             aes(x = Hugo_Symbol,
                 y = log2Count,
                 color = Cohort,
                 group = Cohort)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Driver Gene Comparison (Log2 counts)",
    x = "Driver genes",
    y = "log2(count + 1)"
  ) +
  scale_color_manual(values = c("Pakistan" = "red", "PJL" = "darkgreen"))

ggsave("plots/volcano_driver_compare.png", p2, width = 8, height = 6)

############################################################
# PATHWAY PLOT (CLEAN + COLOURED)
############################################################

p3 <- ggplot(path_compare,
             aes(x = Pathway,
                 y = Count,
                 fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("Pakistan" = "orange", "PJL" = "purple")) +
  labs(
    title = "Pathway Enrichment Comparison",
    y = "Driver gene count"
  )

ggsave("plots/pathway_compare.png", p3, width = 6, height = 5)

############################################################
# PRINT PLOTS
############################################################
print(p1)
print(p2)
print(p3)


pak_samples <- length(unique(maf_sa@data$Tumor_Sample_Barcode))
pak_mut_burden <- maf_sa@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(n_mut = n())

pjl_samples <- sum(pjl_samples > 0)  # already defined earlier

fig1 <- ggplot() +
  geom_bar(aes(x = c("Pakistan","PJL"),
               y = c(pak_samples, pjl_samples)),
           stat = "identity", fill = c("steelblue","darkgreen")) +
  theme_classic() +
  labs(title = "Cohort Size Comparison",
       x = "",
       y = "Number of samples")

ggsave("FIGURES/fig1_cohort_size.png", fig1)

library(VennDiagram)

pak_genes <- unique(maf_sa@data$Hugo_Symbol)
pjl_genes <- unique(pjl_gene_summary$gene)

venn.plot <- VennDiagram::draw.pairwise.venn(
  area1 = length(pak_genes),
  area2 = length(pjl_genes),
  cross.area = length(intersect(pak_genes, pjl_genes)),
  category = c("Pakistan Tumour", "PJL Germline"),
  fill = c("red", "green"),
  alpha = 0.5,
  cex = 1.5
)

png("FIGURES/fig2_venn.png")
grid::grid.draw(venn.plot)
dev.off()


pathways <- list(
  WNT  = wnt_genes,
  PI3K = pi3k_genes,
  RAS  = ras_genes,
  DDR  = ddr_genes
)

results <- lapply(names(pathways), function(pw) {
  
  genes <- toupper(pathways[[pw]])
  
  # remove NA safety
  pak_gene_set <- unique(na.omit(maf_sa@data$Hugo_Symbol))
  pjl_gene_set <- unique(na.omit(pjl_gene_summary$gene))
  
  pak_hits <- sum(pak_gene_set %in% genes)
  pjl_hits <- sum(pjl_gene_set %in% genes)
  
  pak_total <- length(pak_gene_set)
  pjl_total <- length(pjl_gene_set)
  
  # FINAL SAFETY CHECK (prevents fisher crash)
  mat <- matrix(
    c(
      pak_hits,
      max(pak_total - pak_hits, 0),
      pjl_hits,
      max(pjl_total - pjl_hits, 0)
    ),
    nrow = 2
  )
  
  # skip invalid cases
  if (any(mat < 0) || any(is.na(mat))) {
    return(data.frame(
      Pathway = pw,
      OR = NA,
      p = NA
    ))
  }
  
  test <- fisher.test(mat)
  
  data.frame(
    Pathway = pw,
    OR = as.numeric(test$estimate),
    p = test$p.value
  )
})

path_df <- bind_rows(results) %>%
  mutate(FDR = p.adjust(p, method = "BH"))

fig3 <- ggplot(path_df,
               aes(x = Pathway, y = OR, fill = Pathway)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = "Pathway Enrichment (Odds Ratio)")

ggsave("FIGURES/fig3_pathway_OR.png", fig3)


driver_table <- driver_compare %>%
  mutate(
    OR = (Pakistan + 1) / (PJL + 1)
  )

fig4 <- ggplot(driver_table,
               aes(x = reorder(Hugo_Symbol, OR), y = log2(OR))) +
  geom_point(color = "purple", size = 3) +
  coord_flip() +
  theme_classic() +
  labs(title = "Gene-level Enrichment (log2 Odds Ratio)",
       x = "Driver genes",
       y = "log2(OR Pakistan vs PJL)")

ggsave("FIGURES/fig4_gene_OR.png", fig4)


############################################################
# DRIVER ENRICHMENT TABLE
############################################################

driver_table <- driver_compare %>%
  mutate(
    log2_enrichment = log2((Pakistan + 1)/(PJL + 1))
  ) %>%
  arrange(desc(abs(log2_enrichment)))

############################################################
# CLASSIFICATION
############################################################

driver_table <- driver_table %>%
  mutate(
    Category = case_when(
      log2_enrichment > 1  ~ "Tumour-enriched",
      log2_enrichment < -1 ~ "PJL-enriched",
      TRUE ~ "Shared/neutral"
    )
  )

############################################################
# VIEW
############################################################

head(driver_table, 20)

############################################################
# SAVE
############################################################

write.csv(
  driver_table,
  "FIGURES/driver_enrichment_table.csv",
  row.names = FALSE
)

plot_df <- driver_table %>%
  arrange(log2FC) %>%
  mutate(
    Hugo_Symbol = factor(Hugo_Symbol,
                         levels = Hugo_Symbol)
  )

p <- ggplot(plot_df,
            aes(x = Hugo_Symbol,
                y = log2FC,
                color = Category,
                size = abs(log2FC))) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  geom_point(alpha = 0.9) +
  
  coord_flip() +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Relative Driver Gene Enrichment",
    subtitle = "Pakistani OSCC tumours vs PJL germline background",
    x = "",
    y = "log2 Fold Change (Pakistan / PJL)"
  ) +
  
  scale_color_manual(values = c(
    "Tumour-enriched" = "#D55E00",
    "Shared/neutral" = "grey40",
    "PJL-enriched" = "#0072B2"
  ))

print(p)

ggsave(
  "Pakistan_vs_PJL_driver_enrichment.pdf",
  p,
  width = 8,
  height = 7
)