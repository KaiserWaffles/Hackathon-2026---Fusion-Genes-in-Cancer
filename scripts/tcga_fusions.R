#### Load Data
tcgaFusions = read.table('../TCGA/TCGA/cbioportalR_TCGA_fusions.csv',header=TRUE, sep=',')
crc_clinical = read.table('../TCGA/TCGA/clinical/TCGA-COAD_clinical.csv', header = TRUE, sep=',')
panc_clinical = read.table('../TCGA/TCGA/clinical/TCGA-PAAD_clinical.csv', header = TRUE, sep=',')
glio_clinical = read.table('../TCGA/TCGA/clinical/TCGA-GBM_clinical.csv', header = TRUE, sep=',')

### Libraries 
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

# create a new table with only required columns
col_names = c('uniqueSampleKey', 'uniquePatientKey', 'sampleId', 'patientId', 'site1HugoSymbol',"site1EntrezGeneId", 'site2HugoSymbol',"site2EntrezGeneId",'tumorReadCount', 'tumorSplitReadCount', 'tumorPairedEndReadCount','site2EffectOnFrame','eventInfo')
master_fusions = tcgaFusions[,col_names]
master_fusions = master_fusions %>% dplyr::select(patientId, everything())

### Fusion Expressions and Events
table(master_fusions$site2EffectOnFrame, useNA = 'always')


# Summarise read counts by frame status
master_fusions %>%
  filter(!is.na(site2EffectOnFrame)) %>%
  group_by(site2EffectOnFrame) %>%
  summarise(
    n              = n(),
    median_reads   = median(tumorReadCount, na.rm = TRUE),
    median_split   = median(tumorSplitReadCount, na.rm = TRUE),
    median_paired  = median(tumorPairedEndReadCount, na.rm = TRUE)
  ) %>%
  arrange(desc(median_reads))

# using only valid values and better columns, which are tumorReadCount and tumourPairedReadCount
fusions_clean <- master_fusions %>%
  filter(!is.na(site2EffectOnFrame)) %>%
  mutate(
    splitReads  = ifelse(tumorSplitReadCount  < 0, NA, tumorSplitReadCount),
    pairedReads = ifelse(tumorPairedEndReadCount < 0, NA, tumorPairedEndReadCount)
  ) %>%
  filter(!is.na(splitReads) | !is.na(pairedReads))

# pairedRead
pariedRead_effect_boxPlot = fusions_clean %>%
  filter(!is.na(pairedReads)) %>%
  ggplot(aes(x = site2EffectOnFrame, y = pairedReads + 1, fill = site2EffectOnFrame)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  scale_y_log10() +
  labs(title  = 'Paired Read Support by Frame Effect',
       x      = 'Effect on Frame',
       y      = 'Tumour Paired Read Count + 1 (log10)') +
  theme_bw() +
  theme(legend.position = 'none')
pariedRead_effect_boxPlot

####### Fusions and Frequency
# Top 40 fusion pairs by frequency
top_fusions <- master_fusions %>%
  group_by(eventInfo, site2EffectOnFrame) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(eventInfo) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  slice_max(total, n = 40) %>%  
  arrange(desc(total))

# Print to get values for the widget below
print(top_fusions, n = 40)

# Prepare data
fusionFrequency_data <- top_fusions %>%
  mutate(site2EffectOnFrame = ifelse(is.na(site2EffectOnFrame), "Unknown/NA", site2EffectOnFrame))

fusionFrequency_data <- fusionFrequency_data %>%
  mutate(eventInfo = str_remove_all(eventInfo, regex("\\s*fusion", ignore_case = TRUE)))

# remove 'fusion' in eventInfo column
fusionFrequency_data <- fusionFrequency_data %>%
  mutate(eventInfo = str_remove_all(eventInfo, regex("\\s*fusion", ignore_case = TRUE)) %>%
           str_trim())

head(fusionFrequency_data$eventInfo)  # check first five fusions


# categorise effects-on-frame
fusionFrequency_data <- fusionFrequency_data %>%
  mutate(site2EffectOnFrame = factor(site2EffectOnFrame,
                                     levels = c("in-frame", "frameshift", "Unknown/NA")))

fusionFrequency_data <- fusionFrequency_data %>%
  group_by(eventInfo) %>%
  mutate(total = sum(n)) %>%
  ungroup()

fusionFrequency_data <- fusionFrequency_data %>%
  slice_max(total, n = 45)

n_distinct(fusionFrequency_data$eventInfo)  # check how many unique events kept
levels(fusionFrequency_data$site2EffectOnFrame)  # should print all 3 levels




### create a master table for colorectal cancer 
### by merging fuisons and clinical tables 

names(crc_clinical)

# Clinical columns to keep
clinical_cols = c(
  "submitter_id",                    
  "project",
  "tissue_or_organ_of_origin",
  "age_at_diagnosis",
  "primary_diagnosis",
  "ajcc_pathologic_stage",
  "ajcc_pathologic_t",
  "ajcc_pathologic_n",
  "ajcc_pathologic_m",
  "vital_status",
  "days_to_death",
  "days_to_recurrence",
  "progression_or_recurrence",
  "last_known_disease_status",
  "prior_malignancy",
  "prior_treatment",
  "tumor_grade",
  "residual_disease",
  "gender",
  "race",
  "age_at_index",
  "cigarettes_per_day",
  "alcohol_history",
  "follow_ups_disease_response",
  "treatments_pharmaceutical_treatment_type",
  "treatments_pharmaceutical_treatment_outcome",
  "treatments_radiation_treatment_type",
  "treatments_radiation_treatment_outcome"
)

# Subset clinical tables
crc_clinical_subset = crc_clinical[, clinical_cols]
panc_clinical_subset = panc_clinical[, clinical_cols]
glio_clinical_subset = glio_clinical[, clinical_cols]

# Check submitter_id format matches patientId in fusions
head(crc_clinical_subset$submitter_id)
head(master_fusions$patientId)  

# merge fusion and clinical data for colorectal cancer (crc)
master_crc = merge(master_fusions, crc_clinical_subset, by.x = 1, by.y = 1)
master_panc = merge(master_fusions, panc_clinical_subset, by.x = 1, by.y = 1)
master_glio = merge(master_fusions, glio_clinical_subset, by.x = 1, by.y = 1)



### Functions

exp_effect_boxPlot = function(master_table, cancer_type = "Cancer") {
  
  # Step 1 — prepare data
  # Step 1 — prepare data
  fusion_expression <- master_table %>%
    filter(!is.na(site2EffectOnFrame)) %>%
    mutate(
      pairedReads = ifelse(tumorPairedEndReadCount < 0, NA, tumorPairedEndReadCount),
      frameEffect = factor(site2EffectOnFrame, levels = c("in-frame", "frameshift"))
    ) %>%
    filter(!is.na(pairedReads)) %>%
    filter(!is.na(frameEffect)) %>%   
    droplevels()                       
  
  # Step 2 — statistical test (only if both groups present)
  present_levels <- levels(fusion_expression$frameEffect)
  
  if (length(present_levels) == 2) {
    test_result <- wilcox.test(pairedReads ~ frameEffect, data = fusion_expression)
    p_label     <- paste0("Wilcoxon p = ", round(test_result$p.value, 3))
  } else {
    p_label <- paste0("Wilcoxon test skipped — only '", present_levels, "' fusions present")
  }
  
  # Step 3 — summary stats
  stats <- fusion_expression %>%
    group_by(frameEffect) %>%
    summarise(
      n      = n(),
      median = median(pairedReads),
      .groups = "drop"
    )
  print(stats)
  
  # Step 4 — plot
  ggplot(fusion_expression, aes(x = frameEffect, y = pairedReads + 1, fill = frameEffect)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_y_log10() +
    scale_fill_manual(values = c(
      "in-frame"   = "#1D9E75",
      "frameshift" = "#D85A30"
    )) +
    labs(
      title    = paste0("Fusion expression by frame effect — ", cancer_type),
      subtitle = paste0("Tumour paired-end read count as expression proxy | ", p_label),
      x        = NULL,
      y        = "Paired-end read count + 1 (log10)",
      caption  = paste0("n = ", n_distinct(master_table$patientId), " patients")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}

fusionEvent_freq_barChart = function(master_table, top_n = 10, title = NULL) {
  
  # Step 1 — get top N fusion names
  top_events = master_table %>%
    group_by(eventInfo) %>%
    summarise(total = n(), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = top_n) %>%
    pull(eventInfo)
  
  # Step 2 — build frequency table
  fusion_freq_table = master_table %>%
    filter(eventInfo %in% top_events) %>%
    group_by(eventInfo, site2EffectOnFrame) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      site2EffectOnFrame = ifelse(is.na(site2EffectOnFrame), "Unknown/NA", site2EffectOnFrame),
      eventInfo          = str_remove_all(eventInfo, regex("\\s*fusion", ignore_case = TRUE)) %>% str_trim(),
      site2EffectOnFrame = factor(site2EffectOnFrame, levels = c("in-frame", "frameshift", "Unknown/NA"))
    ) %>%
    group_by(eventInfo) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    mutate(eventInfo = fct_reorder(eventInfo, total))
  
  # Step 3 — auto title if not provided
  plot_title <- if (!is.null(title)) title else "Top fusion events by frequency"
  
  # Step 4 — plot
  ggplot(fusion_freq_table, aes(x = n, y = eventInfo, fill = site2EffectOnFrame)) +
    geom_col() +
    scale_fill_manual(values = c(
      "in-frame"   = "#1D9E75",
      "frameshift" = "#D85A30",
      "Unknown/NA" = "#B4B2A9"
    )) +
    labs(
      title    = plot_title,
      subtitle = paste0("n = ", n_distinct(master_table$patientId), " patients"),
      x        = "Number of events",
      y        = NULL,
      fill     = "Frame effect"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position    = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

run_GO_ORA_plot = function(master_table,
                           cancer_type   = "Cancer",
                           frame_filter  = c("in-frame", "frameshift", "both"),
                           gene_partner  = c("both", "site1", "site2"),
                           ont           = "BP",
                           showCategory  = 10,
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.20,
                           universe      = NULL,
                           OrgDb         = org.Hs.eg.db) {
  
  frame_filter = match.arg(frame_filter)
  gene_partner = match.arg(gene_partner)
  
  # Step 1 — filter by frame
  filtered_table = if (frame_filter == "both") master_table else
    master_table %>% filter(site2EffectOnFrame == frame_filter)
  
  # Step 2 — use Entrez IDs directly (already in the data — no bitr needed)
  entrez_ids = switch(gene_partner,
                      "site1" = filtered_table$site1EntrezGeneId,
                      "site2" = filtered_table$site2EntrezGeneId,
                      "both"  = c(filtered_table$site1EntrezGeneId, filtered_table$site2EntrezGeneId)
  ) %>%
    na.omit() %>%
    as.character() %>%
    .[. != ""] %>%
    unique()
  
  partner_label = switch(gene_partner,
                         "site1" = "Gene 1 (5' partner)",
                         "site2" = "Gene 2 (3' partner)",
                         "both"  = "All fusion genes"
  )
  frame_label = switch(frame_filter,
                       "in-frame"   = "In-frame",
                       "frameshift" = "Frameshift",
                       "both"       = "All frames"
  )
  
  cat("Unique Entrez IDs for ORA:", length(entrez_ids), "\n")
  
  if (length(entrez_ids) == 0) {
    warning("No Entrez IDs found.")
    return(NULL)
  }
  
  # Step 3 — ORA (readable = TRUE still works; clusterProfiler maps back to symbols)
  ora_res = enrichGO(
    gene         = entrez_ids,
    universe     = universe,
    OrgDb        = OrgDb,
    keyType      = "ENTREZID",   # explicit now that we skip bitr
    readable     = TRUE,
    ont          = ont,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff
  )
  
  ora_df = as.data.frame(ora_res)
  cat("Enriched GO terms:", nrow(ora_df), "\n")
  print(head(ora_df))
  
  if (nrow(ora_df) == 0) {
    warning("No enriched GO terms. Try ont = 'MF' or 'CC', or relax qvalueCutoff further.")
    return(list(ora_result = ora_res, top_table = NULL,
                dotplot = NULL, cp_barplot = NULL, gg_barplot = NULL))
  }
  
  top_df = head(ora_df[order(ora_df$p.adjust), ], showCategory)
  
  dot_plot = dotplot(ora_res, showCategory = showCategory) +
    labs(
      title    = paste0("GO ", ont, " — ", cancer_type),
      subtitle = paste0(partner_label, " | ", frame_label, " fusions | n = ",
                        n_distinct(master_table$patientId), " patients")
    ) +
    theme_minimal(base_size = 11)
  
  gg_barplot = ggplot(top_df,
                      aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "#D85A30", high = "#B4B2A9", name = "Adjusted\np-value") +
    labs(
      title    = paste0("GO ", ont, " enrichment — ", cancer_type),
      subtitle = paste0(partner_label, " | ", frame_label, " fusions | n = ",
                        n_distinct(master_table$patientId), " patients"),
      x = "Gene count", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())
  
  return(list(
    ora_result = ora_res,
    top_table  = top_df,
    dotplot    = dot_plot,
    cp_barplot = barplot(ora_res, showCategory = showCategory,
                         orderBy = "Count", decreasing = TRUE),
    gg_barplot = gg_barplot
  ))
}

crc_exp_effect_boxPlot = exp_effect_boxPlot(master_crc, cancer_type = "Colorectal Cancer")
crc_exp_effect_boxPlot
panc_exp_effect_boxPlot = exp_effect_boxPlot(master_panc, cancer_type = "Pancreatic Cancer")
panc_exp_effect_boxPlot
glio_exp_effect_boxPlot = exp_effect_boxPlot(master_glio, cancer_type = "Brain Cancer - Glioblastoma")
glio_exp_effect_boxPlot

crc_fusion_freq_chart = fusionEvent_freq_barChart(master_crc, top_n = 10, title = 'Top fusion frequency - Colorectal Cancer')
crc_fusion_freq_chart
panc_fusion_freq_chart = fusionEvent_freq_barChart(master_panc, top_n = 10, title = 'Top fusion frequency - Pancreatic Cancer')
panc_fusion_freq_chart
glio_fusion_freq_chart = fusionEvent_freq_barChart(master_glio, top_n = 10, title = 'Top fusion frequency - Glioblastoma')
glio_fusion_freq_chart

# CRC — site1 (5' partner, typically contributes promoter/regulatory domain)
crc_GO_site1 = run_GO_ORA_plot(master_crc, cancer_type = "Colorectal Cancer",
                               gene_partner = "site1", frame_filter = "both")
crc_GO_site1$dotplot
crc_GO_site1$gg_barplot

# CRC — site2 (3' partner, typically contributes kinase/functional domain)
crc_GO_site2 = run_GO_ORA_plot(master_crc, cancer_type = "Colorectal Cancer",
                               gene_partner = "site2", frame_filter = "both")
crc_GO_site2$dotplot
crc_GO_site2$gg_barplot

# Repeat for other cancers
panc_GO_site1 = run_GO_ORA_plot(master_panc, cancer_type = "Pancreatic Cancer",  gene_partner = "site1")
panc_GO_site1
panc_GO_site2 = run_GO_ORA_plot(master_panc, cancer_type = "Pancreatic Cancer",  gene_partner = "site2")
panc_GO_site2
glio_GO_site1 = run_GO_ORA_plot(master_glio, cancer_type = "Glioblastoma",       gene_partner = "site1")
glio_GO_site1
glio_GO_site2 = run_GO_ORA_plot(master_glio, cancer_type = "Glioblastoma",       gene_partner = "site2")
glio_GO_site2

universe_entrez = c(master_fusions$site1EntrezGeneId,
                    master_fusions$site2EntrezGeneId) %>%
  na.omit() %>% as.character() %>% .[. != ""] %>% unique()

# Check how many unique Entrez IDs are in the CRC gene set
crc_entrez = c(master_crc$site1EntrezGeneId, master_crc$site2EntrezGeneId) %>%
  na.omit() %>% as.character() %>% .[. != ""] %>% unique()

cat("CRC fusion genes:", length(crc_entrez), "\n")
cat("Universe size:", length(universe_entrez), "\n")