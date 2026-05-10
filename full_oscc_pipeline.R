############################################################
# 📦 LIBRARIES
############################################################

library(data.table)
library(dplyr)
library(stringr)
library(maftools)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(uwot)
library(patchwork)
library(TCGAbiolinks)

output_dir <- "FIGURES_new"
dir.create(output_dir, showWarnings = FALSE)

############################################################
# 🎨 PUBLICATION THEME
############################################################

theme_pub <- function(base_size = 14){
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title   = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position   = "right",
      panel.grid.major.y = element_line(color = "grey90"),
      strip.text   = element_text(face = "bold")
    )
}

palette_cohort <- c(
  "Pakistani" = "#D55E00",
  "TCGA_EUR"  = "#0072B2",
  "PJL"       = "#009E73"
)

############################################################
# 🧬 GENE SETS
############################################################

ddr_genes  <- toupper(c("ATM","BRCA1","BRCA2","CHEK2","ATR","CHEK1"))
wnt_genes  <- toupper(c("CTNNB1","APC","AXIN1","AXIN2","TCF7L2","LRP5","LRP6","DKK1"))
pi3k_genes <- toupper(c("PIK3CA","PIK3R1","AKT1","AKT2","PTEN","MTOR","TSC1","TSC2"))
ras_genes  <- toupper(c("KRAS","NRAS","HRAS","BRAF","MAP2K1","MAP2K2","MAPK1"))

driver_genes <- unique(toupper(c(
  "TP53","CDKN2A","NOTCH1","FAT1",
  wnt_genes, pi3k_genes, ras_genes, ddr_genes
)))

noise_genes <- toupper(c(
  "TTN","MUC16","MUC4","MUC5AC","MUC5B",
  "OBSCN","DST","FLG","DNAH1","DNAH2",
  "DNAH5","DNAH8","RYR1","AHNAK2",
  "PCLO","UBR4","LRP1","HMCN2",
  "EPPK1","MUC12","PLEC","SPTBN5",
  "APOB","SBF2","BLTP1","CSMD2"
))

############################################################
# 🧼 HELPER FUNCTIONS
############################################################

clean_gene <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- sub("[,;].*$", "", x)
  x <- trimws(toupper(x))
  x[x %in% c("", ".", "NA", "N/A")] <- NA
  return(x)
}

safe_col <- function(df, col, default = NA){
  if(col %in% colnames(df)) df[[col]] else rep(default, nrow(df))
}

############################################################
# 🔄 ANNOVAR → MAF
############################################################

convert_annovar_to_maf <- function(file){
  
  df <- fread(file, data.table = FALSE)
  
  required <- c("Chr","Start","End","Ref","Alt",
                "Gene.ensGene","ExonicFunc.ensGene")
  stopifnot(all(required %in% colnames(df)))
  
  df$Chr   <- as.character(df$Chr)
  df$Start <- as.numeric(df$Start)
  df$End   <- as.numeric(df$End)
  df$Ref   <- as.character(df$Ref)
  df$Alt   <- as.character(df$Alt)
  
  # drop invalid rows
  df <- df[
    !is.na(df$Gene.ensGene) & df$Gene.ensGene != "" &
      !is.na(df$Chr) & !is.na(df$Start) &
      !is.na(df$Ref) & !is.na(df$Alt), ]
  
  gene <- clean_gene(df$Gene.ensGene)
  
  maf <- data.frame(
    Hugo_Symbol          = gene,
    Chromosome           = gsub("^chr", "", df$Chr),
    Start_Position       = df$Start,
    End_Position         = df$End,
    Reference_Allele     = toupper(df$Ref),
    Tumor_Seq_Allele2    = toupper(df$Alt),
    Variant_Classification = df$ExonicFunc.ensGene,
    Variant_Type         = ifelse(nchar(df$Ref)==1 & nchar(df$Alt)==1, "SNP", "INDEL"),
    Tumor_Sample_Barcode = tools::file_path_sans_ext(basename(file)),
    HGVSp_Short          = sapply(strsplit(df$AAChange.ensGene, ":"),
                                  function(x) tail(x, 1)),
    CADD     = suppressWarnings(as.numeric(safe_col(df, "CADD_phred"))),
    SIFT     = safe_col(df, "SIFT_pred"),
    PolyPhen = safe_col(df, "Polyphen2_HDIV_pred"),
    SAS_AF   = suppressWarnings(as.numeric(safe_col(df, "1000g2015aug_sas"))),
    ALL_AF   = suppressWarnings(as.numeric(safe_col(df, "1000g2015aug_all"))),
    Cytoband = safe_col(df, "cytoBand"),
    COSMIC   = safe_col(df, "cosmic70"),
    Func     = safe_col(df, "Func.ensGene"),
    stringsAsFactors = FALSE
  )
  
  # drop rows with no valid gene
  maf <- maf[!is.na(maf$Hugo_Symbol) &
               maf$Hugo_Symbol != "NA" &
               maf$Hugo_Symbol != "", ]
  
  # standardise variant classifications
  maf$Variant_Classification <- dplyr::case_when(
    maf$Variant_Classification %in%
      c("nonsynonymous SNV","nonsynonymous_SNV")      ~ "Missense_Mutation",
    maf$Variant_Classification == "stopgain"           ~ "Nonsense_Mutation",
    maf$Variant_Classification == "stoploss"           ~ "Nonstop_Mutation",
    maf$Variant_Classification %in%
      c("frameshift deletion","frameshift_deletion")   ~ "Frame_Shift_Del",
    maf$Variant_Classification %in%
      c("frameshift insertion","frameshift_insertion") ~ "Frame_Shift_Ins",
    maf$Variant_Classification %in%
      c("nonframeshift deletion","nonframeshift_deletion") ~ "In_Frame_Del",
    maf$Variant_Classification %in%
      c("nonframeshift insertion","nonframeshift_insertion") ~ "In_Frame_Ins",
    maf$Variant_Classification == "splicing"           ~ "Splice_Site",
    TRUE ~ NA_character_
  )
  
  maf <- maf[!is.na(maf$Variant_Classification), ]
  return(maf)
}

############################################################
# 📁 LOAD + MERGE ANNOVAR FILES
############################################################

files <- list.files("annovar_detailed",
                    pattern = "\\.hg38_multianno\\.txt$",
                    full.names = TRUE)

maf_df <- dplyr::bind_rows(lapply(files, convert_annovar_to_maf)) %>%
  dplyr::filter(!is.na(Hugo_Symbol),
                Hugo_Symbol != "NA",
                Hugo_Symbol != "")

############################################################
# 📊 BUILD MAF OBJECTS
# NOTE: three clearly named objects, no overwriting
#   maf_pk_raw  = read.maf of everything
#   maf_pk_full = noise filtered (used for TMB + lollipops)
#   maf_pk      = noise + driver filtered (used for comparisons)
############################################################

maf_pk_raw  <- read.maf(maf = maf_df)

maf_pk_full <- subsetMaf(
  maf   = maf_pk_raw,
  query = paste0("!Hugo_Symbol %in% c('",
                 paste(noise_genes, collapse = "','"), "')")
)

maf_pk <- subsetMaf(
  maf   = maf_pk_full,   # ✅ subset from full, not from itself
  genes = driver_genes
)

# add cohort labels
maf_pk_full@data[, Cohort := "Pakistani"]
maf_pk@data[,      Cohort := "Pakistani"]

cat("Pakistani samples (full):", 
    length(unique(maf_pk_full@data$Tumor_Sample_Barcode)), "\n")
cat("Pakistani samples (driver):", 
    length(unique(maf_pk@data$Tumor_Sample_Barcode)), "\n")

print(getGeneSummary(maf_pk_full)[1:15, ])

############################################################
# 🧬 TCGA — LOAD, FILTER, CLEAN
############################################################

tcga <- tcgaLoad("HNSC", source = "MC3")

oral_sites <- c(
  "Oral Tongue","Base of tongue","Floor of mouth",
  "Buccal Mucosa","Alveolar Ridge","Hard Palate",
  "Oral Cavity","Lip"
)

keep_samples <- tcga@clinical.data %>%
  dplyr::filter(
    tolower(race) == "white",
    anatomic_neoplasm_subdivision %in% oral_sites
  ) %>%
  dplyr::pull(Tumor_Sample_Barcode)

tcga_eur <- subsetMaf(tcga, tsb = keep_samples)

# clean TCGA MAF
tcga_eur@data <- data.table::as.data.table(tcga_eur@data)
tcga_eur@data <- tcga_eur@data[, !duplicated(names(tcga_eur@data)), with = FALSE]
tcga_eur@data[, Hugo_Symbol := toupper(as.character(Hugo_Symbol))]
tcga_eur@data <- tcga_eur@data[
  !is.na(Hugo_Symbol) & Hugo_Symbol != "" & Hugo_Symbol != "NA"
]

valid_classes <- c(
  "Missense_Mutation","Nonsense_Mutation",
  "Frame_Shift_Del","Frame_Shift_Ins",
  "In_Frame_Del","In_Frame_Ins",
  "Splice_Site","Nonstop_Mutation"
)
tcga_eur@data <- tcga_eur@data[Variant_Classification %in% valid_classes]

tcga_eur_full <- subsetMaf(
  maf   = tcga_eur,
  query = paste0("!Hugo_Symbol %in% c('",
                 paste(noise_genes, collapse = "','"), "')")
)
tcga_eur_full@data[, Cohort := "TCGA_EUR"]

tcga_eur <- subsetMaf(tcga_eur_full, genes = driver_genes)
tcga_eur@data[, Cohort := "TCGA_EUR"]

cat("\nTCGA samples (full):",
    length(unique(tcga_eur_full@data$Tumor_Sample_Barcode)), "\n")
print(getGeneSummary(tcga_eur)[1:15, ])

############################################################
# 🔗 MERGE
############################################################

maf_combined_df <- data.table::rbindlist(
  list(maf_pk@data, tcga_eur@data),
  fill = TRUE
)

maf_combined <- maftools::read.maf(maf = maf_combined_df)

cohort_map <- unique(maf_combined@data[, .(Tumor_Sample_Barcode, Cohort)])

cat("\n=== COHORT DISTRIBUTION ===\n")
print(table(maf_combined@data$Cohort))

############################################################
# 📊 TMB — FIXED
# Uses maf_pk_full and tcga_eur_full (NOT driver-filtered)
# total_perMB is already in mut/Mb — DO NOT multiply by 1e6
############################################################

capture_size <- 38

tmb_pk   <- tmb(maf_pk_full,   captureSize = capture_size)
tmb_pk$Cohort <- "Pakistani"

tmb_tcga <- tmb(tcga_eur_full, captureSize = capture_size)
tmb_tcga$Cohort <- "TCGA_EUR"

tmb_combined <- rbind(tmb_pk, tmb_tcga)

# QC check — these should be sensible numbers
cat("\nPakistani TMB summary:\n");   print(summary(tmb_pk$total_perMB))
cat("\nTCGA TMB summary:\n");        print(summary(tmb_tcga$total_perMB))

median_df <- tmb_combined %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarise(MedianTMB = median(total_perMB, na.rm = TRUE))

p_tmb <- ggplot(tmb_combined,
                aes(x = Cohort, y = total_perMB, fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.4, width = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  geom_point(data = median_df,
             aes(x = Cohort, y = MedianTMB),
             shape = 23, size = 3, fill = "white") +
  geom_text(data = median_df,
            aes(x = Cohort, y = MedianTMB,
                label = round(MedianTMB, 2)),
            vjust = -1.2, size = 4) +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(
    method  = "wilcox.test",
    label   = "p.format",   # shows exact value
    label.x = 1.5,
    label.y = 47
  ) +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "Tumor Mutation Burden Comparison",
       y = "Mutations per Mb", x = NULL)

ggsave(file.path(output_dir, "Fig_TMB_Final.pdf"),
       p_tmb, width = 6, height = 5)

############################################################
# 📊 MUTATION SPECTRUM
############################################################

mut_prop <- maf_combined@data %>%
  dplyr::filter(!is.na(Variant_Classification)) %>%
  dplyr::group_by(Cohort, Variant_Classification) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  dplyr::group_by(Cohort) %>%
  dplyr::mutate(Proportion = Count / sum(Count))

mutation_colors <- c(
  "Missense_Mutation" = "#4DBBD5",
  "Nonsense_Mutation" = "#E64B35",
  "Frame_Shift_Del"   = "#00A087",
  "Frame_Shift_Ins"   = "#3C5488",
  "In_Frame_Del"      = "#F39B7F",
  "In_Frame_Ins"      = "#8491B4",
  "Splice_Site"       = "#91D1C2",
  "Nonstop_Mutation"  = "#DC0000"
)

p_mut <- ggplot(mut_prop,
                aes(x = Cohort, y = Proportion,
                    fill = Variant_Classification)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = mutation_colors) +
  theme_pub() +
  labs(title = "Mutation Spectrum by Cohort",
       y = "Proportion of Mutations", fill = "Variant Type")

ggsave(file.path(output_dir, "Fig_MutationSpectrum.pdf"),
       p_mut, width = 6, height = 5)

############################################################
# 📊 DRIVER GENE PANEL
############################################################

driver_table <- maf_combined@data %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::group_by(Hugo_Symbol, Cohort) %>%
  dplyr::summarise(MutatedSamples = n_distinct(Tumor_Sample_Barcode),
                   .groups = "drop") %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::filter(max(MutatedSamples) >= 2) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!Hugo_Symbol %in%
                  c("AKT1","NRAS","KRAS","DKK1","BRAF","CHEK1","TCF7L2"))

gene_order <- driver_table %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(total = sum(MutatedSamples)) %>%
  dplyr::arrange(total) %>%
  dplyr::pull(Hugo_Symbol)

driver_table$Hugo_Symbol <- factor(driver_table$Hugo_Symbol, levels = gene_order)

p_driver <- ggplot(driver_table,
                   aes(x = Hugo_Symbol, y = MutatedSamples, fill = Cohort)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.75), width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = palette_cohort) +
  theme_pub(base_size = 15) +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 11),
        plot.margin = margin(10, 20, 10, 10)) +
  labs(title = "Driver Gene Mutations",
       y = "Mutated Samples", x = NULL)

# multi-panel with TMB
p_multi <- ggarrange(p_tmb, p_driver,
                     ncol = 2, widths = c(1, 2),
                     common.legend = TRUE, legend = "top")

ggsave(file.path(output_dir, "Fig_TMB_Driver_Panel.pdf"),
       p_multi, width = 12, height = 8)

############################################################
# 📊 ONCOPLOTS
############################################################

top_genes <- getGeneSummary(maf_combined) %>%
  dplyr::arrange(desc(MutatedSamples)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Hugo_Symbol)

pdf(file.path(output_dir, "Fig_Oncoplot_Pakistani.pdf"), 10, 7)
oncoplot(maf_pk, genes = top_genes, removeNonMutated = FALSE)
dev.off()

pdf(file.path(output_dir, "Fig_Oncoplot_TCGA.pdf"), 12, 7)
oncoplot(tcga_eur, genes = top_genes, removeNonMutated = FALSE)
dev.off()

############################################################
# 📊 FISHER TEST + VOLCANO + FOREST
############################################################

total_pk   <- length(unique(maf_pk@data$Tumor_Sample_Barcode))
total_tcga <- length(unique(tcga_eur@data$Tumor_Sample_Barcode))

gene_table <- maf_combined@data %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::group_by(Hugo_Symbol, Cohort) %>%
  dplyr::summarise(n = n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Cohort, values_from = n, values_fill = 0) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    fisher_out = list(fisher.test(matrix(
      c(Pakistani, total_pk   - Pakistani,
        TCGA_EUR,  total_tcga - TCGA_EUR), nrow = 2))),
    OR      = fisher_out$estimate,
    p.value = fisher_out$p.value
  ) %>%
  dplyr::select(-fisher_out) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(p.value)

write.table(gene_table,
            file.path(output_dir, "Table_Fisher.tsv"),
            sep = "\t", row.names = FALSE)

plot_df <- gene_table %>%
  dplyr::filter(!is.na(OR), OR > 0, Pakistani >= 2 | TCGA_EUR >= 2) %>%
  dplyr::mutate(log2OR = log2(OR), neglog10p = -log10(p.value))

label_df <- plot_df %>%
  dplyr::filter(FDR < 0.05 | abs(log2OR) > 1) %>%
  dplyr::arrange(p.value) %>%
  dplyr::slice_head(n = 15)

p_volcano <- ggplot(plot_df, aes(x = log2OR, y = neglog10p)) +
  geom_point(aes(color = FDR < 0.05), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey70","red")) +
  geom_text_repel(data = label_df, aes(label = Hugo_Symbol),
                  size = 3, max.overlaps = 100, box.padding = 0.3) +
  theme_pub() +
  labs(title = "Driver Gene Enrichment",
       x = "log2 Odds Ratio", y = "-log10(p-value)") +
  theme(legend.position = "none")

forest_df <- plot_df %>%
  dplyr::arrange(FDR) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    SE      = sqrt(1/Pakistani + 1/(total_pk   - Pakistani) +
                     1/TCGA_EUR  + 1/(total_tcga - TCGA_EUR)),
    CI_low  = exp(log(OR) - 1.96*SE),
    CI_high = exp(log(OR) + 1.96*SE)
  )

p_forest <- ggplot(forest_df, aes(reorder(Hugo_Symbol, OR), OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high)) +
  coord_flip() + scale_y_log10() +
  theme_pub() + labs(title = "Forest Plot")

ggsave(file.path(output_dir, "Fig_Volcano_Forest.pdf"),
       p_volcano + p_forest, width = 12, height = 5)

############################################################
# 📊 JACCARD
############################################################

n_pk   <- length(unique(maf_pk@data$Tumor_Sample_Barcode))
n_tcga <- length(unique(tcga_eur@data$Tumor_Sample_Barcode))

freq_pk <- maf_pk@data %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(pct = n_distinct(Tumor_Sample_Barcode) / n_pk * 100)

freq_tcga <- tcga_eur@data %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(pct = n_distinct(Tumor_Sample_Barcode) / n_tcga * 100)

genes_pk_5    <- freq_pk   %>% dplyr::filter(pct >= 5) %>% dplyr::pull(Hugo_Symbol)
genes_tcga_5  <- freq_tcga %>% dplyr::filter(pct >= 5) %>% dplyr::pull(Hugo_Symbol)
top20_pk      <- freq_pk   %>% dplyr::arrange(desc(pct)) %>% dplyr::slice_head(n=20) %>% dplyr::pull(Hugo_Symbol)
top20_tcga    <- freq_tcga %>% dplyr::arrange(desc(pct)) %>% dplyr::slice_head(n=20) %>% dplyr::pull(Hugo_Symbol)

jaccard_5pct  <- length(intersect(genes_pk_5,  genes_tcga_5))  / length(union(genes_pk_5,  genes_tcga_5))
jaccard_top20 <- length(intersect(top20_pk,    top20_tcga))    / length(union(top20_pk,    top20_tcga))
jaccard_driver <- length(intersect(
  c("TP53","FAT1","CDKN2A","NOTCH1","PIK3CA"),
  c("TP53","FAT1","CDKN2A","NOTCH1","CASP8","PIK3CA"))) /
  length(union(
    c("TP53","FAT1","CDKN2A","NOTCH1","PIK3CA"),
    c("TP53","FAT1","CDKN2A","NOTCH1","CASP8","PIK3CA")))

cat("Jaccard (>=5%):", jaccard_5pct, "\n")
cat("Jaccard (top20):", jaccard_top20, "\n")
cat("Jaccard (drivers):", jaccard_driver, "\n")

############################################################
# 📊 LOLLIPOP PLOTS
############################################################

lolli_dir <- file.path(output_dir, "Lollipop_Plots")
dir.create(lolli_dir, showWarnings = FALSE)

make_lollipop <- function(maf_obj, gene, width = 12, height = 5){
  
  pdf(file.path(lolli_dir, paste0(gene, "_lollipop.pdf")),
      width = width, height = height)
  par(mar = c(5,4,5,2), xpd = TRUE)
  
  tryCatch({
    lollipopPlot(
      maf            = maf_obj,
      gene           = gene,
      AACol          = "HGVSp_Short",
      showMutationRate = TRUE,
      showDomainLabel  = TRUE,
      labelPos       = "all",
      repel          = TRUE,
      labPosAngle    = 60,
      labPosSize     = 0.7
    )
    title(main = paste0(gene, " — Pakistani OSCC Cohort"),
          font.main = 2, cex.main = 1.4, line = 2)
  }, error = function(e){
    message("Skipping ", gene, ": ", e$message)
  })
  
  dev.off()
}

for(gene in c("TP53","FAT1","NOTCH1","MTOR")){
  make_lollipop(maf_pk_full, gene)
}

############################################################
# 📊 UMAP
############################################################

mut_mat <- maf_combined@data %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  dplyr::mutate(value = 1) %>%
  tidyr::pivot_wider(names_from  = Hugo_Symbol,
                     values_from = value,
                     values_fill = list(value = 0)) %>%
  dplyr::arrange(Tumor_Sample_Barcode)

sample_ids  <- mut_mat$Tumor_Sample_Barcode
mut_mat_num <- mut_mat %>% dplyr::select(-Tumor_Sample_Barcode)
mut_mat_num[is.na(mut_mat_num)] <- 0
mut_mat_num <- data.frame(lapply(mut_mat_num, as.numeric))

set.seed(42)
umap_res <- uwot::umap(mut_mat_num, n_neighbors = 10,
                       min_dist = 0.3, metric = "hamming")

umap_df <- data.frame(UMAP1 = umap_res[,1],
                      UMAP2 = umap_res[,2],
                      Sample = sample_ids) %>%
  merge(cohort_map, by.x = "Sample", by.y = "Tumor_Sample_Barcode")

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Cohort)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.5) +
  scale_color_manual(values = palette_cohort) +
  theme_pub() + labs(title = "UMAP — Driver Gene Mutational Profiles")

ggsave(file.path(output_dir, "Fig_UMAP.pdf"), p_umap, width = 6.5, height = 5.5)

############################################################
# 📊 HIGH-CONFIDENCE FUNCTIONAL VARIANTS
############################################################

high_conf <- maf_combined@data %>%
  dplyr::filter(
    Func %in% c("exonic","splicing"),
    Variant_Classification != "Silent",
    !is.na(CADD), CADD >= 20,
    SIFT %in% c("D"),
    PolyPhen %in% c("D","P"),
    is.na(SAS_AF) | SAS_AF < 0.01
  ) %>%
  dplyr::arrange(desc(CADD))

write.table(high_conf,
            file.path(output_dir, "Table_HighConfidenceVariants.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# CADD distribution plot
p_cadd <- ggplot(high_conf, aes(x = CADD, fill = Cohort)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "Distribution of Predicted Deleteriousness",
       x = "CADD Phred Score", y = "Variant Count")

ggsave(file.path(output_dir, "Fig_CADD_Distribution.pdf"),
       p_cadd, width = 7, height = 5)

############################################################
# 📁 TABLE EXPORTS
############################################################

write.table(getGeneSummary(maf_combined),
            file.path(output_dir, "Table_TopGenes.tsv"),
            sep = "\t", row.names = FALSE)

write.table(driver_table,
            file.path(output_dir, "Table_DriverGenes.tsv"),
            sep = "\t", row.names = FALSE)

write.table(getSampleSummary(maf_combined) %>%
              merge(cohort_map, by = "Tumor_Sample_Barcode"),
            file.path(output_dir, "Table_TMB.tsv"),
            sep = "\t", row.names = FALSE)

cat("\nDone. All figures saved to:", output_dir, "\n")


# Check raw mutation counts before TMB calculation
cat("=== RAW MUTATION COUNTS ===\n")
cat("Pakistani total mutations:", nrow(maf_pk_full@data), "\n")
cat("Pakistani samples:", length(unique(maf_pk_full@data$Tumor_Sample_Barcode)), "\n")
cat("Pakistani mutations per sample (raw):\n")
print(summary(table(maf_pk_full@data$Tumor_Sample_Barcode)))

cat("\nTCGA total mutations:", nrow(tcga_eur_full@data), "\n")
cat("TCGA samples:", length(unique(tcga_eur_full@data$Tumor_Sample_Barcode)), "\n")
cat("TCGA mutations per sample (raw):\n")
print(summary(table(tcga_eur_full@data$Tumor_Sample_Barcode)))

# Check what tmb() is actually returning before any scaling
cat("\n=== RAW TMB() OUTPUT ===\n")
cat("Pakistani total_perMB range:\n")
print(summary(tmb_pk$total_perMB))
cat("\nTCGA total_perMB range:\n")
print(summary(tmb_tcga$total_perMB))

# Check column names in tmb output
cat("\nTMB column names:\n")
print(colnames(tmb_pk))



# Figure 1 — overall TMB (all genes, noise filtered only)
tmb_pk_all   <- tmb(maf_pk_full,   captureSize = 38)
tmb_tcga_all <- tmb(tcga_eur_full, captureSize = 38)

# Figure 2A — driver gene TMB
tmb_pk_driver   <- tmb(maf_pk,   captureSize = 38)
tmb_tcga_driver <- tmb(tcga_eur, captureSize = 38)

# Verify all four medians
cat("Overall Pakistani median:",     median(tmb_pk_all$total_perMB), "\n")
cat("Overall TCGA median:",          median(tmb_tcga_all$total_perMB), "\n")
cat("Driver Pakistani median:",      median(tmb_pk_driver$total_perMB), "\n")
cat("Driver TCGA median:",           median(tmb_tcga_driver$total_perMB), "\n")



# Driver gene mutation burden — mutations per sample (not per Mb)
driver_burden_pk <- maf_pk@data %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::summarise(driver_muts = n_distinct(Hugo_Symbol)) %>%
  dplyr::mutate(Cohort = "Pakistani")

driver_burden_tcga <- tcga_eur@data %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::summarise(driver_muts = n_distinct(Hugo_Symbol)) %>%
  dplyr::mutate(Cohort = "TCGA_EUR")

driver_burden <- rbind(driver_burden_pk, driver_burden_tcga)

# Check medians — these should be close to 5 and 2
cat("Driver Pakistani median:", 
    median(driver_burden_pk$driver_muts), "\n")
cat("Driver TCGA median:",      
    median(driver_burden_tcga$driver_muts), "\n")

# Wilcoxon test
wt <- wilcox.test(
  driver_burden_pk$driver_muts,
  driver_burden_tcga$driver_muts
)
cat("Driver burden p-value:", 
    format(wt$p.value, scientific = TRUE), "\n")


############################################################
# 📊 FIGURE 2: TMB MULTI-PANEL
# Panel A: Overall TMB (violin)
# Panel B: Driver gene mutation burden (boxplot)
# Panel C: Top mutated gene comparison (bar)
############################################################

# --- Panel A: Overall TMB ---

tmb_all <- rbind(
  tmb_pk_all   %>% dplyr::mutate(Cohort = "Pakistani"),
  tmb_tcga_all %>% dplyr::mutate(Cohort = "TCGA_EUR")
)

median_all <- tmb_all %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarise(MedianTMB = median(total_perMB, na.rm = TRUE))

p_tmb_overall <- ggplot(tmb_all,
                        aes(x = Cohort, y = total_perMB, fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.4, width = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  geom_point(data = median_all,
             aes(x = Cohort, y = MedianTMB),
             shape = 23, size = 3, fill = "white") +
  geom_text(data = median_all,
            aes(x = Cohort, y = MedianTMB,
                label = round(MedianTMB, 2)),
            vjust = -1.2, size = 4) +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(method = "wilcox.test",
                     label    = "p.format",
                     label.x  = 1.5,
                     label.y  = 47) +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "A: Tumor Mutation Burden",
       y = "Mutations per Mb", x = NULL)

# --- Panel B: Driver Gene Mutation Burden ---

median_driver <- driver_burden %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarise(MedianDriver = median(driver_muts, na.rm = TRUE))

p_tmb_driver <- ggplot(driver_burden,
                       aes(x = Cohort, y = driver_muts, fill = Cohort)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  geom_point(data = median_driver,
             aes(x = Cohort, y = MedianDriver),
             shape = 23, size = 3, fill = "white") +
  geom_text(data = median_driver,
            aes(x = Cohort, y = MedianDriver,
                label = round(MedianDriver, 1)),
            vjust = -1.2, size = 4) +
  stat_compare_means(method = "wilcox.test",
                     label   = "p.format",
                     label.x = 1.5,
                     label.y = max(driver_burden$driver_muts) * 0.95) +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "B: Driver Gene Mutations per Tumour",
       y = "Mutated Driver Genes per Sample", x = NULL)

# --- Panel C: Top Mutated Gene Comparison ---

# Top genes by total mutations across both cohorts
top_genes_combined <- maf_combined@data %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(total = n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Hugo_Symbol)

driver_freq <- maf_combined@data %>%
  dplyr::filter(Hugo_Symbol %in% top_genes_combined) %>%
  dplyr::group_by(Hugo_Symbol, Cohort) %>%
  dplyr::summarise(MutatedSamples = n_distinct(Tumor_Sample_Barcode),
                   .groups = "drop")

# order genes by Pakistani frequency for readability
gene_order <- driver_freq %>%
  dplyr::filter(Cohort == "Pakistani") %>%
  dplyr::arrange(MutatedSamples) %>%
  dplyr::pull(Hugo_Symbol)

driver_freq$Hugo_Symbol <- factor(driver_freq$Hugo_Symbol,
                                  levels = gene_order)

p_genes <- ggplot(driver_freq,
                  aes(x = Hugo_Symbol,
                      y = MutatedSamples,
                      fill = Cohort)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.75),
           width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = palette_cohort) +
  theme_pub(base_size = 13) +
  theme(
    axis.text.y  = element_text(size = 11, face = "bold.italic"),
    axis.text.x  = element_text(size = 10),
    plot.margin  = margin(10, 20, 10, 10)
  ) +
  labs(title = "C: Top Mutated Driver Genes",
       y = "Mutated Samples", x = NULL)

# --- Combine all three panels ---

p_fig2 <- ggarrange(
  ggarrange(p_tmb_overall, p_tmb_driver,
            ncol = 2, nrow = 1,
            widths = c(1, 1),
            common.legend = TRUE,
            legend = "none"),
  p_genes,
  nrow = 2,
  heights = c(1, 1.4),
  common.legend = TRUE,
  legend = "top"
)

ggsave(
  file.path(output_dir, "Fig2_TMB_Driver_Genes_Panel.pdf"),
  p_fig2,
  width  = 12,
  height = 14
)