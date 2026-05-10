############################################################
# 📦 LIBRARIES
############################################################
library(data.table)
library(dplyr)
library(maftools)

############################################################
# 🧬 GENE SETS
############################################################
ddr_genes <- toupper(c("ATM","BRCA1","BRCA2","CHEK2","ATR","CHEK1"))
wnt_genes <- toupper(c("CTNNB1","APC","AXIN1","AXIN2","TCF7L2","LRP5","LRP6","DKK1"))
pi3k_genes <- toupper(c("PIK3CA","PIK3R1","AKT1","AKT2","PTEN","MTOR","TSC1","TSC2"))
ras_genes <- toupper(c("KRAS","NRAS","HRAS","BRAF","MAP2K1","MAP2K2","MAPK1"))

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
# 🧼 GENE CLEANING
############################################################
clean_gene <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- sub("[,;].*$", "", x)
  x <- trimws(x)
  x <- toupper(x)
  x[x %in% c("", ".", "NA", "N/A")] <- NA
  return(x)
}

############################################################
# 🔄 ANNOVAR → MAF (STANDARDIZED)
############################################################
convert_annovar_to_maf <- function(file){
  
  df <- fread(file, data.table = FALSE)
  
  required <- c("Chr","Start","End","Ref","Alt","Gene.ensGene","ExonicFunc.ensGene")
  stopifnot(all(required %in% colnames(df)))
  
  # Clean gene
  gene <- clean_gene(df$Gene.ensGene)
  
  maf <- data.frame(
    Hugo_Symbol = gene,
    Chromosome = gsub("^chr", "", df$Chr),
    Start_Position = as.numeric(df$Start),
    End_Position = as.numeric(df$End),
    Reference_Allele = toupper(df$Ref),
    Tumor_Seq_Allele2 = toupper(df$Alt),
    Variant_Classification = df$ExonicFunc.ensGene,
    Variant_Type = ifelse(nchar(df$Ref)==1 & nchar(df$Alt)==1, "SNP", "INDEL"),
    Tumor_Sample_Barcode = tools::file_path_sans_ext(basename(file)),
    stringsAsFactors = FALSE
  )
  
  # Clean rows
  maf <- maf %>%
    filter(!is.na(Hugo_Symbol))
  
  # Standardize classification
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
    filter(!is.na(Variant_Classification))
  
  return(maf)
}

############################################################
# 📁 LOAD BOTH COHORTS (IDENTICAL PROCESSING)
############################################################

process_cohort <- function(path){
  files <- list.files(path,
                      pattern="\\.hg38_multianno\\.txt$",
                      full.names=TRUE)
  
  maf_list <- lapply(files, convert_annovar_to_maf)
  maf_df <- bind_rows(maf_list)
  
  maf_df <- maf_df %>%
    filter(!is.na(Hugo_Symbol))
  
  return(maf_df)
}

############################################################
# 🧬 BUILD DATASETS
############################################################

maf_pak_raw <- process_cohort("annovar_results/pakistan")
maf_ind_raw <- process_cohort("annovar_results/india")

############################################################
# 🧹 REMOVE NOISE (OPTIONAL BUT CONSISTENT)
############################################################

maf_pak_clean <- maf_pak_raw %>%
  filter(!Hugo_Symbol %in% noise_genes)

maf_ind_clean <- maf_ind_raw %>%
  filter(!Hugo_Symbol %in% noise_genes)

############################################################
# 🧬 NON-SYNONYMOUS FILTER (IDENTICAL)
############################################################

nonsynonymous <- c(
  "Missense_Mutation","Nonsense_Mutation",
  "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site"
)

maf_pak_clean <- maf_pak_clean %>%
  filter(Variant_Classification %in% nonsynonymous)

maf_ind_clean <- maf_ind_clean %>%
  filter(Variant_Classification %in% nonsynonymous)

############################################################
# 🧬 GENE SETS FOR JACCARD
############################################################

pak_genes_all <- unique(maf_pak_clean$Hugo_Symbol)
ind_genes_all <- unique(maf_ind_clean$Hugo_Symbol)

############################################################
# 🔁 RECURRENT GENES ≥5%
############################################################

get_recurrent <- function(maf_df){
  maf_df %>%
    group_by(Hugo_Symbol) %>%
    summarise(freq = n_distinct(Tumor_Sample_Barcode)) %>%
    mutate(freq_pct = freq / length(unique(maf_df$Tumor_Sample_Barcode)) * 100) %>%
    filter(freq_pct >= 5) %>%
    pull(Hugo_Symbol)
}

pak_recurrent <- get_recurrent(maf_pak_clean)
ind_recurrent <- get_recurrent(maf_ind_clean)

############################################################
# 🔝 TOP 20 GENES
############################################################

get_top20 <- function(maf_df){
  maf_df %>%
    count(Hugo_Symbol, sort = TRUE) %>%
    slice_head(n = 20) %>%
    pull(Hugo_Symbol)
}

pak_top20 <- get_top20(maf_pak_clean)
ind_top20 <- get_top20(maf_ind_clean)

############################################################
# 🧬 DRIVER GENES
############################################################

pak_driver <- intersect(pak_genes_all, driver_genes)
ind_driver <- intersect(ind_genes_all, driver_genes)

pak_recurrent <- get_recurrent(maf_pak_clean)
ind_recurrent <- get_recurrent(maf_ind_clean)

length(pak_recurrent)
length(ind_recurrent)

############################################################
# 📊 JACCARD FUNCTION
############################################################

jaccard_index <- function(a, b){
  length(intersect(a,b)) / length(union(a,b))
}

############################################################
# 📈 FINAL RESULTS
############################################################

jaccard_results <- data.frame(
  Category = c("All genes","Recurrent ≥5%","Top 20","Driver"),
  Jaccard = c(
    jaccard_index(pak_genes_all, ind_genes_all),
    jaccard_index(pak_recurrent, ind_recurrent),
    jaccard_index(pak_top20, ind_top20),
    jaccard_index(pak_driver, ind_driver)
  )
)

print(jaccard_results)