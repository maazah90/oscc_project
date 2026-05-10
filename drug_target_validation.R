############################################################
# 🧬 CLINICAL TARGET VALIDATION (FIXED + UPGRADED)
############################################################

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

############################################################
# 📌 SAFE GENE CLEANING (CRITICAL FIX)
############################################################
clean_gene <- function(x){
  x <- as.character(x)
  x <- sub("[,;].*$", "", x)   # keep first gene only
  x <- toupper(trimws(x))
  x[x %in% c("", ".", "NA")] <- NA
  return(x)
}

############################################################
# 🔄 ANNOVAR → CLEAN TABLE
############################################################
convert_annovar_to_clean <- function(file){
  
  df <- fread(file)
  
  df$Gene <- clean_gene(df$Gene.ensGene)
  df$Func <- df$ExonicFunc.ensGene
  
  df <- df %>%
    filter(!is.na(Gene)) %>%
    filter(Func != "synonymous SNV")
  
  sample <- gsub("\\.hg38_multianno.*", "", basename(file))
  
  data.frame(
    Hugo_Symbol = df$Gene,
    Variant_Class = df$Func,
    Sample = sample,
    stringsAsFactors = FALSE
  )
}

files <- list.files("annovar_results",
                    pattern = "\\.hg38_multianno\\.txt$",
                    full.names = TRUE)

maf_df <- bind_rows(lapply(files, convert_annovar_to_clean))

############################################################
# 🧬 VARIANT STANDARDISATION
############################################################
maf_df$Variant_Classification <- case_when(
  maf_df$Variant_Class %in% c("nonsynonymous SNV","nonsynonymous_SNV") ~ "Missense_Mutation",
  maf_df$Variant_Class == "stopgain" ~ "Nonsense_Mutation",
  maf_df$Variant_Class == "stoploss" ~ "Nonstop_Mutation",
  maf_df$Variant_Class %in% c("frameshift deletion","frameshift_deletion") ~ "Frame_Shift_Del",
  maf_df$Variant_Class %in% c("frameshift insertion","frameshift_insertion") ~ "Frame_Shift_Ins",
  maf_df$Variant_Class %in% c("nonframeshift deletion","nonframeshift_deletion") ~ "In_Frame_Del",
  maf_df$Variant_Class %in% c("nonframeshift insertion","nonframeshift_insertion") ~ "In_Frame_Ins",
  TRUE ~ NA_character_
)

maf_df <- maf_df %>% filter(!is.na(Variant_Classification))


############################################################
# 🧪 DRUG CLASS MAPPING
############################################################

drug_targets <- list(
  EGFR_family = c("EGFR","ERBB2","ERBB3"),
  PI3K_pathway = c("PIK3CA","PIK3R1","PTEN","MTOR"),
  MAPK_pathway = c("KRAS","NRAS","HRAS","BRAF","MAP2K1"),
  Cell_cycle = c("CDK4","CDK6","CCND1","CDKN2A"),
  DNA_damage = c("BRCA1","BRCA2","ATM","ATR"),
  RTK = c("FGFR1","FGFR2","FGFR3","ALK","ROS1","RET","MET")
)

maf_df$Drug_Target_Class <- NA

for (i in names(drug_targets)) {
  maf_df$Drug_Target_Class[
    maf_df$Hugo_Symbol %in% drug_targets[[i]]
  ] <- i
}

############################################################
# 📊 TABLE 1: ACTIONABLE GENES
############################################################
maf_df <- maf_df %>%
  mutate(Actionable = Hugo_Symbol %in% actionable_genes)

clinical_targets <- maf_df %>%
  filter(Actionable) %>%
  group_by(Hugo_Symbol, Drug_Target_Class) %>%
  summarise(
    Mutated_Samples = n_distinct(Sample),
    Total_Mutations = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Mutated_Samples))

write.table(
  clinical_targets,
  "Table_Clinical_Actionable.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 📊 PATIENT-LEVEL ACTIONABILITY (FIXED — CORRECT)
############################################################

# 1. Subset actionable mutations
actionable_df <- maf_combined@data %>%
  filter(Hugo_Symbol %in% actionable_genes)

# 2. Patients with actionable mutations
patients_actionable <- actionable_df %>%
  distinct(Tumor_Sample_Barcode, Cohort)

# 3. Total patients
total_patients <- maf_combined@data %>%
  distinct(Tumor_Sample_Barcode, Cohort)

# 4. Create matching keys (SAFE version)
patients_actionable_key <- patients_actionable %>%
  mutate(key = paste(Tumor_Sample_Barcode, Cohort))

total_patients <- total_patients %>%
  mutate(key = paste(Tumor_Sample_Barcode, Cohort))

# 5. Compute actionable rate
actionable_rate <- total_patients %>%
  mutate(Has_Actionable = key %in% patients_actionable_key$key) %>%
  group_by(Cohort) %>%
  summarise(
    Total = n(),
    Actionable = sum(Has_Actionable),
    Percent = Actionable / Total * 100,
    .groups = "drop"
  )
############################################################
# 📊 FIGURE: ACTIONABLE GENE FREQUENCY
############################################################

p_actionable <- ggplot(clinical_targets,
                       aes(x = reorder(Hugo_Symbol, Mutated_Samples),
                           y = Mutated_Samples,
                           fill = Drug_Target_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Clinically Actionable Mutations",
    x = NULL,
    y = "Mutated Samples"
  )

ggsave("Fig_Actionable_Genes.pdf", p_actionable, width = 7, height = 6)

############################################################
# 📊 OPTIONAL: NETWORK (CLEANED)
############################################################

library(igraph)
library(ggraph)

edges <- clinical_targets %>%
  filter(!is.na(Drug_Target_Class)) %>%
  select(Hugo_Symbol, Drug_Target_Class)

g <- graph_from_data_frame(edges)

pdf("Fig_Drug_Gene_Network.pdf", 8, 6)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.5) +
  geom_node_point(size = 4) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

dev.off()

############################################################
# 🧬 CLINICAL TARGET VALIDATION: COHORT COMPARISON
############################################################

############################################################
# 🎯 ACTIONABLE GENE SET (same as before)
############################################################
actionable_genes <- toupper(c(
  "EGFR","ERBB2","ERBB3",
  "PIK3CA","PIK3R1","PTEN","MTOR",
  "KRAS","NRAS","HRAS","BRAF","MAP2K1",
  "CDK4","CDK6","CCND1","CDKN2A",
  "BRCA1","BRCA2","ATM","ATR",
  "FGFR1","FGFR2","FGFR3",
  "ALK","ROS1","RET","MET"
))

############################################################
# 🧬 SUBSET ACTIONABLE FROM COMBINED DATA
############################################################

actionable_df <- maf_combined@data %>%
  filter(Hugo_Symbol %in% actionable_genes)

############################################################
# 📊 TABLE 1: ACTIONABLE GENES PER COHORT
############################################################

actionable_summary <- actionable_df %>%
  group_by(Hugo_Symbol, Cohort) %>%
  summarise(
    Mutated_Samples = n_distinct(Tumor_Sample_Barcode),
    .groups = "drop"
  )

write.table(
  actionable_summary,
  file.path(output_dir, "Table_Actionable_ByCohort.tsv"),
  sep = "\t",
  row.names = FALSE
)

############################################################
# 📊 PATIENT-LEVEL ACTIONABILITY (CRITICAL)
############################################################

# Patients with ≥1 actionable mutation
patients_actionable <- actionable_df %>%
  distinct(Tumor_Sample_Barcode, Cohort)

# Total patients
total_patients <- maf_combined@data %>%
  distinct(Tumor_Sample_Barcode, Cohort)

actionable_rate <- total_patients %>%
  mutate(Has_Actionable = Tumor_Sample_Barcode %in% patients_actionable$Tumor_Sample_Barcode) %>%
  group_by(Cohort) %>%
  summarise(
    Total = n(),
    Actionable = sum(Has_Actionable),
    Percent = Actionable / Total * 100,
    .groups = "drop"
  )

write.table(
  actionable_rate,
  file.path(output_dir, "Table_Actionable_Rate_ByCohort.tsv"),
  sep = "\t",
  row.names = FALSE
)

############################################################
# 📊 STATISTICAL TEST: IS ACTIONABILITY DIFFERENT?
############################################################

# contingency table
sa_total <- actionable_rate$Total[actionable_rate$Cohort=="South_Asian"]
tcga_total <- actionable_rate$Total[actionable_rate$Cohort=="TCGA_EUR"]

sa_act <- actionable_rate$Actionable[actionable_rate$Cohort=="South_Asian"]
tcga_act <- actionable_rate$Actionable[actionable_rate$Cohort=="TCGA_EUR"]

fisher_res <- fisher.test(matrix(c(
  sa_act,
  sa_total - sa_act,
  tcga_act,
  tcga_total - tcga_act
), nrow = 2))

write.table(
  data.frame(
    OddsRatio = fisher_res$estimate,
    P_value = fisher_res$p.value
  ),
  file.path(output_dir, "Table_Actionable_Fisher.tsv"),
  sep = "\t",
  row.names = FALSE
)

# AFTER fisher_res is created



############################################################
# 📊 FIGURE 1: ACTIONABLE PATIENT PERCENTAGE
############################################################
p_actionable_rate <- ggplot(actionable_rate,
                            aes(x = Cohort,
                                y = Percent,
                                fill = Cohort)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(
    title = "Patients with Actionable Mutations",
    y = "Percentage (%)"
  )

p_actionable_rate <- p_actionable_rate +
  annotate("text",
           x = mean(as.numeric(factor(actionable_rate$Cohort))),
           y = max(actionable_rate$Percent) * 1.05,
           label = paste0("p = ", signif(fisher_res$p.value, 3)),
           size = 5)

ggsave(file.path(output_dir, "Fig_Actionable_Rate.pdf"),
       p_actionable_rate, width = 5, height = 5)


############################################################
# 📊 FIGURE 2: ACTIONABLE GENE COMPARISON
############################################################

p_actionable_genes <- ggplot(actionable_summary,
                             aes(x = reorder(Hugo_Symbol, Mutated_Samples),
                                 y = Mutated_Samples,
                                 fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(
    title = "Clinically Actionable Genes (SA vs TCGA)",
    y = "Mutated Samples"
  )

ggsave(file.path(output_dir, "Fig_Actionable_Genes_Comparison.pdf"),
       p_actionable_genes, width = 7, height = 6)

############################################################
# 📊 FIGURE 3: PATHWAY-LEVEL ACTIONABILITY
############################################################

drug_targets <- list(
  EGFR_family = c("EGFR","ERBB2","ERBB3"),
  PI3K_pathway = c("PIK3CA","PIK3R1","PTEN","MTOR"),
  MAPK_pathway = c("KRAS","NRAS","HRAS","BRAF","MAP2K1"),
  Cell_cycle = c("CDK4","CDK6","CCND1","CDKN2A"),
  DNA_damage = c("BRCA1","BRCA2","ATM","ATR"),
  RTK = c("FGFR1","FGFR2","FGFR3","ALK","ROS1","RET","MET")
)

actionable_df$Pathway <- NA

for (i in names(drug_targets)) {
  actionable_df$Pathway[
    actionable_df$Hugo_Symbol %in% drug_targets[[i]]
  ] <- i
}

pathway_actionable <- actionable_df %>%
  filter(!is.na(Pathway)) %>%
  group_by(Cohort, Pathway) %>%
  summarise(
    Mutated_Samples = n_distinct(Tumor_Sample_Barcode),
    .groups = "drop"
  )

p_pathway_actionable <- ggplot(pathway_actionable,
                               aes(x = Pathway,
                                   y = Mutated_Samples,
                                   fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = palette_cohort) +
  theme_pub() +
  labs(
    title = "Actionable Pathway Burden",
    y = "Mutated Samples"
  )

ggsave(file.path(output_dir, "Fig_Actionable_Pathways.pdf"),
       p_pathway_actionable, width = 7, height = 5)
