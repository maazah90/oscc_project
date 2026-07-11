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

output_dir <- "FIGURES"
dir.create(output_dir, showWarnings = FALSE)

############################################################
# PUBLICATION THEME
############################################################

theme_pub <- function(base_size = 14){
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      legend.title       = element_text(face = "bold"),
      legend.position    = "right",
      panel.grid.major.y = element_line(color = "grey90"),
      strip.text         = element_text(face = "bold")
    )
}

# ✅ Indian cohort added to palette
palette_cohort <- c(
  "Pakistani" = "#D55E00",   # orange (your original)
  "TCGA_EUR"  = "#0072B2",   # blue (your original)
  "Indian"    = "#009E73",   # teal-green (your original)
  "PJL"       = "#CC79A7"    # pink/mauve — distinct from all above
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
# 🔧 HELPER FUNCTIONS
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

# Recurrent genes in ≥ pct% of samples (fixed scoping)
get_recurrent <- function(maf_df, pct = 5){
  n_samples <- dplyr::n_distinct(maf_df$Tumor_Sample_Barcode)
  maf_df %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(freq = dplyr::n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
    dplyr::mutate(freq_pct = freq / n_samples * 100) %>%
    dplyr::filter(freq_pct >= pct) %>%
    dplyr::pull(Hugo_Symbol)
}

get_top20 <- function(maf_df){
  maf_df %>%
    dplyr::count(Hugo_Symbol, sort = TRUE) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::pull(Hugo_Symbol)
}

jaccard_index <- function(a, b){
  length(intersect(a, b)) / length(union(a, b))
}

############################################################
# 🔄 ANNOVAR → MAF
# Single shared function used for BOTH Pakistani and Indian cohorts.
# Pulls all extra columns (CADD, SIFT, PolyPhen, etc.) needed
# downstream for high-confidence variant filtering and lollipops.
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
  
  df <- df[
    !is.na(df$Gene.ensGene) & df$Gene.ensGene != "" &
      !is.na(df$Chr) & !is.na(df$Start) &
      !is.na(df$Ref) & !is.na(df$Alt), ]
  
  gene <- clean_gene(df$Gene.ensGene)
  
  maf <- data.frame(
    Hugo_Symbol            = gene,
    Chromosome             = gsub("^chr", "", df$Chr),
    Start_Position         = df$Start,
    End_Position           = df$End,
    Reference_Allele       = toupper(df$Ref),
    Tumor_Seq_Allele2      = toupper(df$Alt),
    Variant_Classification = df$ExonicFunc.ensGene,
    Variant_Type           = ifelse(nchar(df$Ref)==1 & nchar(df$Alt)==1, "SNP", "INDEL"),
    Tumor_Sample_Barcode   = tools::file_path_sans_ext(basename(file)),
    HGVSp_Short            = sapply(strsplit(df$AAChange.ensGene, ":"),
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
  
  maf <- maf[!is.na(maf$Hugo_Symbol) &
               maf$Hugo_Symbol != "NA" &
               maf$Hugo_Symbol != "", ]
  
  maf$Variant_Classification <- dplyr::case_when(
    maf$Variant_Classification %in%
      c("nonsynonymous SNV","nonsynonymous_SNV")          ~ "Missense_Mutation",
    maf$Variant_Classification == "stopgain"               ~ "Nonsense_Mutation",
    maf$Variant_Classification == "stoploss"               ~ "Nonstop_Mutation",
    maf$Variant_Classification %in%
      c("frameshift deletion","frameshift_deletion")       ~ "Frame_Shift_Del",
    maf$Variant_Classification %in%
      c("frameshift insertion","frameshift_insertion")     ~ "Frame_Shift_Ins",
    maf$Variant_Classification %in%
      c("nonframeshift deletion","nonframeshift_deletion") ~ "In_Frame_Del",
    maf$Variant_Classification %in%
      c("nonframeshift insertion","nonframeshift_insertion") ~ "In_Frame_Ins",
    maf$Variant_Classification == "splicing"               ~ "Splice_Site",
    TRUE ~ NA_character_
  )
  
  maf <- maf[!is.na(maf$Variant_Classification), ]
  return(maf)
}

############################################################
# 📁 LOAD ANNOVAR FILES — SHARED LOADER
############################################################

load_annovar_cohort <- function(path){
  files <- list.files(path,
                      pattern    = "\\.hg38_multianno\\.txt$",
                      full.names = TRUE)
  if(length(files) == 0) stop(paste("No ANNOVAR files found in:", path))
  dplyr::bind_rows(lapply(files, convert_annovar_to_maf)) %>%
    dplyr::filter(!is.na(Hugo_Symbol),
                  Hugo_Symbol != "NA",
                  Hugo_Symbol != "")
}

############################################################
# 🇮🇳 INDIAN COHORT — LOAD FIRST
# Loaded before Pakistani so objects are available for merging.
# Mirrors the exact same three-tier structure as Pakistani cohort:
#   maf_ind_raw  = read.maf of everything
#   maf_ind_full = noise filtered (used for TMB + lollipops)
#   maf_ind      = noise + driver filtered (used for comparisons)
############################################################

maf_df_indian <- load_annovar_cohort("annovar_indian_cohort")

maf_ind_raw  <- maftools::read.maf(maf = maf_df_indian, isTCGA = FALSE, verbose = FALSE)

maf_ind_full <- subsetMaf(
  maf   = maf_ind_raw,
  query = paste0("!Hugo_Symbol %in% c('",
                 paste(noise_genes, collapse = "','"), "')")
)

maf_ind <- subsetMaf(
  maf   = maf_ind_full,
  genes = driver_genes
)

maf_ind_full@data[, Cohort := "Indian"]
maf_ind@data[,      Cohort := "Indian"]

cat("Indian samples (full):  ", length(unique(maf_ind_full@data$Tumor_Sample_Barcode)), "\n")
cat("Indian samples (driver):", length(unique(maf_ind@data$Tumor_Sample_Barcode)), "\n")
print(getGeneSummary(maf_ind_full)[1:15, ])

############################################################
# 🇵🇰 PAKISTANI COHORT
############################################################

maf_df_pakistani <- load_annovar_cohort("annovar_detailed")

maf_pk_raw  <- maftools::read.maf(maf = maf_df_pakistani)

maf_pk_full <- subsetMaf(
  maf   = maf_pk_raw,
  query = paste0("!Hugo_Symbol %in% c('",
                 paste(noise_genes, collapse = "','"), "')")
)

maf_pk <- subsetMaf(
  maf   = maf_pk_full,
  genes = driver_genes
)

maf_pk_full@data[, Cohort := "Pakistani"]
maf_pk@data[,      Cohort := "Pakistani"]

cat("Pakistani samples (full):  ", length(unique(maf_pk_full@data$Tumor_Sample_Barcode)), "\n")
cat("Pakistani samples (driver):", length(unique(maf_pk@data$Tumor_Sample_Barcode)), "\n")
print(getGeneSummary(maf_pk_full)[1:15, ])

############################################################
# 🧬 TCGA — LOAD, FILTER, CLEAN
############################################################

options(timeout = 360)
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

cat("\nTCGA samples (full):  ", length(unique(tcga_eur_full@data$Tumor_Sample_Barcode)), "\n")
cat("TCGA samples (driver):", length(unique(tcga_eur@data$Tumor_Sample_Barcode)), "\n")
print(getGeneSummary(tcga_eur)[1:15, ])

############################################################
# 🔀 MERGE ALL THREE COHORTS
############################################################

maf_combined_df <- data.table::rbindlist(
  list(maf_pk@data, maf_ind@data, tcga_eur@data),
  fill = TRUE
)

maf_combined <- maftools::read.maf(maf = maf_combined_df)

cohort_map <- unique(maf_combined@data[, .(Tumor_Sample_Barcode, Cohort)])

cat("\n=== COHORT DISTRIBUTION ===\n")
print(table(maf_combined@data$Cohort))

############################################################
# 📊 TMB
# Uses *_full objects (noise filtered only, NOT driver filtered)
############################################################

capture_size <- 38

tmb_pk   <- tmb(maf_pk_full,   captureSize = capture_size); tmb_pk$Cohort   <- "Pakistani"
tmb_ind  <- tmb(maf_ind_full,  captureSize = capture_size); tmb_ind$Cohort  <- "Indian"
tmb_tcga <- tmb(tcga_eur_full, captureSize = capture_size); tmb_tcga$Cohort <- "TCGA_EUR"

tmb_combined <- rbind(tmb_pk, tmb_ind, tmb_tcga)

cat("\nPakistani TMB summary:\n"); print(summary(tmb_pk$total_perMB))
cat("\nIndian TMB summary:\n");    print(summary(tmb_ind$total_perMB))
cat("\nTCGA TMB summary:\n");      print(summary(tmb_tcga$total_perMB))

median_df <- tmb_combined %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarise(MedianTMB = median(total_perMB, na.rm = TRUE))

p_tmb <- ggplot(tmb_combined, aes(x = Cohort, y = total_perMB, fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.4, width = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  geom_point(data = median_df, aes(x = Cohort, y = MedianTMB),
             shape = 23, size = 3, fill = "white") +
  geom_text(data = median_df, aes(x = Cohort, y = MedianTMB,
                                  label = round(MedianTMB, 2)), vjust = -1.2, size = 4) +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(method = "kruskal.test",   # 3-group → Kruskal-Wallis
                     label  = "p.format",
                     label.x = 1.5, label.y = 47) +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "Tumor Mutation Burden Comparison",
       y = "Mutations per Mb", x = NULL)

ggsave(file.path(output_dir, "Fig_TMB_Final_3_cohorts.pdf"), p_tmb, width = 7, height = 5)

############################################################
# 🔬 PAIRWISE TMB (post-hoc Wilcoxon)
############################################################

pairwise.wilcox.test(
  tmb_combined$total_perMB,
  tmb_combined$Cohort,
  p.adjust.method = "BH"
)

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

p_mut <- ggplot(mut_prop, aes(x = Cohort, y = Proportion,
                              fill = Variant_Classification)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = mutation_colors) +
  theme_pub() +
  labs(title = "Mutation Spectrum by Cohort",
       y = "Proportion of Mutations", fill = "Variant Type")

ggsave(file.path(output_dir, "Fig_MutationSpectrum_3_cohort.pdf"), p_mut, width = 7, height = 5)

############################################################
# 🧬 DRIVER GENE PANEL
############################################################

driver_table <- maf_combined@data %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::group_by(Hugo_Symbol, Cohort) %>%
  dplyr::summarise(MutatedSamples = n_distinct(Tumor_Sample_Barcode), .groups = "drop") %>%
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
  labs(title = "Driver Gene Mutations", y = "Mutated Samples", x = NULL)

p_multi <- ggarrange(p_tmb, p_driver, ncol = 2, widths = c(1, 2),
                     common.legend = TRUE, legend = "top")
ggsave(file.path(output_dir, "Fig_TMB_Driver_Panel_3_cohort.pdf"), p_multi, width = 14, height = 8)

############################################################
# 🗺️ ONCOPLOTS
############################################################

top_genes <- getGeneSummary(maf_combined) %>%
  dplyr::arrange(desc(MutatedSamples)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(Hugo_Symbol)

pdf(file.path(output_dir, "Fig_Oncoplot_Pakistani.pdf"), 10, 7)
oncoplot(maf_pk, genes = top_genes, removeNonMutated = FALSE)
dev.off()

pdf(file.path(output_dir, "Fig_Oncoplot_Indian.pdf"), 10, 7)
oncoplot(maf_ind, genes = top_genes, removeNonMutated = FALSE)
dev.off()

pdf(file.path(output_dir, "Fig_Oncoplot_TCGA.pdf"), 12, 7)
oncoplot(tcga_eur, genes = top_genes, removeNonMutated = FALSE)
dev.off()

############################################################
# 🔬 FISHER TEST — ALL PAIRWISE
############################################################

total_pk   <- length(unique(maf_pk@data$Tumor_Sample_Barcode))
total_ind  <- length(unique(maf_ind@data$Tumor_Sample_Barcode))
total_tcga <- length(unique(tcga_eur@data$Tumor_Sample_Barcode))

run_fisher_pair <- function(data_a, data_b, total_a, total_b,
                            label_a, label_b){
  combined <- dplyr::bind_rows(
    data_a %>% dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
      dplyr::group_by(Hugo_Symbol) %>%
      dplyr::summarise(n = n_distinct(Tumor_Sample_Barcode), .groups="drop") %>%
      dplyr::mutate(Cohort = label_a),
    data_b %>% dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
      dplyr::group_by(Hugo_Symbol) %>%
      dplyr::summarise(n = n_distinct(Tumor_Sample_Barcode), .groups="drop") %>%
      dplyr::mutate(Cohort = label_b)
  ) %>%
    tidyr::pivot_wider(names_from = Cohort, values_from = n, values_fill = 0) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      fisher_out = list(fisher.test(matrix(
        c(.data[[label_a]], total_a - .data[[label_a]],
          .data[[label_b]], total_b - .data[[label_b]]), nrow = 2))),
      OR      = fisher_out$estimate,
      p.value = fisher_out$p.value
    ) %>%
    dplyr::select(-fisher_out) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH"),
                  Comparison = paste(label_a, "vs", label_b)) %>%
    dplyr::arrange(p.value)
  combined
}

fisher_pk_ind  <- run_fisher_pair(as.data.frame(maf_pk@data),
                                  as.data.frame(maf_ind@data),
                                  total_pk, total_ind,
                                  "Pakistani", "Indian")

fisher_pk_tcga <- run_fisher_pair(as.data.frame(maf_pk@data),
                                  as.data.frame(tcga_eur@data),
                                  total_pk, total_tcga,
                                  "Pakistani", "TCGA_EUR")

fisher_ind_tcga <- run_fisher_pair(as.data.frame(maf_ind@data),
                                   as.data.frame(tcga_eur@data),
                                   total_ind, total_tcga,
                                   "Indian", "TCGA_EUR")

fisher_all <- dplyr::bind_rows(fisher_pk_ind, fisher_pk_tcga, fisher_ind_tcga)

write.table(fisher_all, file.path(output_dir, "Table_Fisher_3way.tsv"),
            sep = "\t", row.names = FALSE)

############################################################
# 🌋 VOLCANO (Pakistani vs TCGA — primary comparison)
############################################################

plot_df <- fisher_pk_tcga %>%
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
  labs(title = "Driver Gene Enrichment: Pakistani vs TCGA",
       x = "log2 Odds Ratio", y = "-log10(p-value)") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "Fig_Volcano_PK_TCGA.pdf"), p_volcano, width = 7, height = 5)

############################################################
# 📐 JACCARD — THREE-WAY
############################################################

# Use *_full data frames for all three (not driver-subsetted)
pk_df   <- as.data.frame(maf_pk_full@data)
ind_df  <- as.data.frame(maf_ind_full@data)
tcga_df <- as.data.frame(tcga_eur_full@data)

pk_genes_all   <- unique(pk_df$Hugo_Symbol)
ind_genes_all  <- unique(ind_df$Hugo_Symbol)
tcga_genes_all <- unique(tcga_df$Hugo_Symbol)

pk_recurrent   <- get_recurrent(pk_df)
ind_recurrent  <- get_recurrent(ind_df)
tcga_recurrent <- get_recurrent(tcga_df)

pk_top20   <- get_top20(pk_df)
ind_top20  <- get_top20(ind_df)
tcga_top20 <- get_top20(tcga_df)

pk_driver   <- intersect(pk_genes_all,   driver_genes)
ind_driver  <- intersect(ind_genes_all,  driver_genes)
tcga_driver <- intersect(tcga_genes_all, driver_genes)

jaccard_3way <- data.frame(
  Comparison = c("Pakistan vs India", "Pakistan vs TCGA", "India vs TCGA"),
  All_Genes = c(
    jaccard_index(pk_genes_all,  ind_genes_all),
    jaccard_index(pk_genes_all,  tcga_genes_all),
    jaccard_index(ind_genes_all, tcga_genes_all)
  ),
  Recurrent_5pct = c(
    jaccard_index(pk_recurrent,  ind_recurrent),
    jaccard_index(pk_recurrent,  tcga_recurrent),
    jaccard_index(ind_recurrent, tcga_recurrent)
  ),
  Top20 = c(
    jaccard_index(pk_top20,  ind_top20),
    jaccard_index(pk_top20,  tcga_top20),
    jaccard_index(ind_top20, tcga_top20)
  ),
  Driver = c(
    jaccard_index(pk_driver,  ind_driver),
    jaccard_index(pk_driver,  tcga_driver),
    jaccard_index(ind_driver, tcga_driver)
  )
)

cat("\n=== THREE-WAY JACCARD INDICES ===\n")
print(jaccard_3way)
write.table(jaccard_3way, file.path(output_dir, "Table_Jaccard_3way.tsv"),
            sep = "\t", row.names = FALSE)

############################################################
# 🎞️ LOLLIPOP PLOTS (Pakistani + Indian)
############################################################

lolli_dir <- file.path(output_dir, "Lollipop_Plots")
dir.create(lolli_dir, showWarnings = FALSE)

make_lollipop <- function(maf_obj, gene, cohort_label, width = 12, height = 5){
  pdf(file.path(lolli_dir, paste0(gene, "_", cohort_label, "_lollipop.pdf")),
      width = width, height = height)
  par(mar = c(5,4,5,2), xpd = TRUE)
  tryCatch({
    lollipopPlot(maf = maf_obj, gene = gene, AACol = "HGVSp_Short",
                 showMutationRate = TRUE, showDomainLabel = TRUE,
                 labelPos = "all", repel = TRUE,
                 labPosAngle = 60, labPosSize = 0.7)
    title(main = paste0(gene, " — ", cohort_label, " OSCC Cohort"),
          font.main = 2, cex.main = 1.4, line = 2)
  }, error = function(e) message("Skipping ", gene, ": ", e$message))
  dev.off()
}

for(gene in c("TP53","FAT1","NOTCH1","MTOR")){
  make_lollipop(maf_pk_full,  gene, "Pakistani")
  make_lollipop(maf_ind_full, gene, "Indian")
}

############################################################
# 🗺️ UMAP — ALL THREE COHORTS
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

umap_df <- data.frame(UMAP1 = umap_res[,1], UMAP2 = umap_res[,2],
                      Sample = sample_ids) %>%
  merge(cohort_map, by.x = "Sample", by.y = "Tumor_Sample_Barcode")

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Cohort)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.5) +
  scale_color_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "UMAP — Driver Gene Mutational Profiles (3 Cohorts)")

ggsave(file.path(output_dir, "Fig_UMAP_3_cohort.pdf"), p_umap, width = 7, height = 5.5)

############################################################
# 🏆 HIGH-CONFIDENCE FUNCTIONAL VARIANTS
############################################################

high_conf <- maf_combined@data %>%
  dplyr::filter(
    Func %in% c("exonic","splicing"),
    Variant_Classification != "Silent",
    !is.na(CADD), CADD >= 20,
    SIFT     %in% c("D"),
    PolyPhen %in% c("D","P"),
    is.na(SAS_AF) | SAS_AF < 0.01
  ) %>%
  dplyr::arrange(desc(CADD))

write.table(high_conf,
            file.path(output_dir, "Table_HighConfidenceVariants.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p_cadd <- ggplot(high_conf, aes(x = CADD, fill = Cohort)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(title = "Distribution of Predicted Deleteriousness",
       x = "CADD Phred Score", y = "Variant Count")

ggsave(file.path(output_dir, "Fig_CADD_Distribution.pdf"), p_cadd, width = 7, height = 5)

############################################################
# 📋 TABLE EXPORTS
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



############################################################
# 📊 DRIVER GENE FREQUENCY CORRELATION — Pakistani vs Indian
############################################################

# Calculate mutation frequency (%) per driver gene in each cohort
n_pk  <- dplyr::n_distinct(maf_pk_full@data$Tumor_Sample_Barcode)
n_ind <- dplyr::n_distinct(maf_ind_full@data$Tumor_Sample_Barcode)

freq_pk_driver <- as.data.frame(maf_pk_full@data) %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(freq_pk = n_distinct(Tumor_Sample_Barcode) / n_pk * 100,
                   .groups = "drop")

freq_ind_driver <- as.data.frame(maf_ind_full@data) %>%
  dplyr::filter(Hugo_Symbol %in% driver_genes) %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(freq_ind = n_distinct(Tumor_Sample_Barcode) / n_ind * 100,
                   .groups = "drop")

# Merge — keep all driver genes, fill 0 if absent in either cohort
freq_corr <- dplyr::full_join(freq_pk_driver, freq_ind_driver,
                              by = "Hugo_Symbol") %>%
  tidyr::replace_na(list(freq_pk = 0, freq_ind = 0))

# Correlation statistics
cor_test <- cor.test(freq_corr$freq_pk, freq_corr$freq_ind,
                     method = "spearman")  # Spearman — robust to small n and non-normality

r_val <- round(cor_test$estimate, 3)
p_val <- signif(cor_test$p.value, 3)

cat("Spearman r:", r_val, "\n")
cat("p-value:   ", p_val, "\n")

# Label only genes mutated in ≥10% of either cohort
label_genes <- freq_corr %>%
  dplyr::filter(freq_pk >= 10 | freq_ind >= 10)

# Plot
p_corr <- ggplot(freq_corr, aes(x = freq_pk, y = freq_ind)) +
  
  # diagonal reference line (perfect agreement)
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey60", linewidth = 0.5) +
  
  # regression line with CI
  geom_smooth(method = "lm", se = TRUE,
              color = "#D55E00", fill = "#D55E0033", linewidth = 0.8) +
  
  # points — size scaled by average frequency
  geom_point(aes(size = (freq_pk + freq_ind) / 2),
             color = "#333333", alpha = 0.75) +
  
  # label high-frequency genes only
  ggrepel::geom_text_repel(
    data        = label_genes,
    aes(label   = Hugo_Symbol),
    size        = 3.5,
    fontface    = "bold.italic",
    box.padding = 0.4,
    max.overlaps = 20
  ) +
  
  # correlation annotation
  annotate("text",
           x     = max(freq_corr$freq_pk) * 0.05,
           y     = max(freq_corr$freq_ind) * 0.95,
           label = paste0("Spearman r = ", r_val,
                          "\np = ", p_val),
           hjust = 0, size = 4.5, fontface = "bold") +
  
  scale_size_continuous(range = c(2, 6), guide = "none") +
  
  theme_pub() +
  labs(
    title    = "Driver Gene Mutation Frequency Correlation",
    subtitle = "Pakistani vs Indian OSCC Cohorts",
    x        = "Pakistani Cohort — Mutation Frequency (%)",
    y        = "Indian Cohort — Mutation Frequency (%)"
  )

ggsave(file.path(output_dir, "Fig_FreqCorrelation_PK_IND.pdf"),
       p_corr, width = 7, height = 6)