suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(patchwork)
  library(scales)
})


OUTDIR <- "fusion_plots"
dir.create(OUTDIR, showWarnings = FALSE)

save_plot <- function(p, filename, w = 12, h = 8) {
  ggsave(file.path(OUTDIR, filename), p,
         width = w, height = h, dpi = 150, bg = "white")
  message("  Saved → ", filename)
}


save_km <- function(p_km, filename, w = 2800, h = 2500) {
  png(file.path(OUTDIR, filename), width = w, height = h, res = 200, bg = "white")
  print(p_km)
  dev.off()
  message("  Saved → ", filename)
}


theme_tcga <- theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    strip.background = element_rect(fill = "grey93"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )
theme_set(theme_tcga)

CANCER_COLORS <- c(
  BRCA     = "#E41A1C", SARC     = "#377EB8", COADREAD = "#4DAF4A",
  LUAD     = "#FF7F00", LUSC     = "#984EA3", OV       = "#A65628",
  BLCA     = "#F781BF", SKCM     = "#999999", STAD     = "#66C2A5",
  LIHC     = "#FC8D62", HNSC     = "#8DA0CB", LGG      = "#E78AC3",
  UCEC     = "#A6D854", GBM      = "#E5C494", ESCA     = "#FFD92F",
  CESC     = "#B3B3FF", UCS      = "#FF99CC", KIRC     = "#CCFF99",
  ACC      = "#FFD700", PRAD     = "#00CED1"
)


message("\n── 1. Loading data ──")

fusions  <- read_csv("cbioportalR_TCGA_fusions.csv",            show_col_types = FALSE)
clinical <- read_csv("clinical/clinical_shared_columns.csv",    show_col_types = FALSE)

brca_clin     <- read_csv("clinical/TCGA-BRCA_clinical.csv",   show_col_types = FALSE)
sarc_clin     <- read_csv("clinical/TCGA-SARC_clinical.csv",   show_col_types = FALSE)
coad_clin     <- read_csv("clinical/TCGA-COAD_clinical.csv",   show_col_types = FALSE)
read_clin_df  <- read_csv("clinical/TCGA-READ_clinical.csv",   show_col_types = FALSE)
luad_clin     <- read_csv("clinical/TCGA-LUAD_clinical.csv",   show_col_types = FALSE)

# COAD + READ are merged as "coadread" in the fusion file
# Coerce all columns to character first to avoid type conflicts between the two files
coadread_clin <- bind_rows(
  coad_clin    %>% mutate(across(everything(), as.character)),
  read_clin_df %>% mutate(across(everything(), as.character))
)



# studyId → short cancer abbreviation
study_map <- tribble(
  ~studyId,                               ~cancer_type,
  "brca_tcga_pan_can_atlas_2018",         "BRCA",
  "sarc_tcga_pan_can_atlas_2018",         "SARC",
  "coadread_tcga_pan_can_atlas_2018",     "COADREAD",
  "luad_tcga_pan_can_atlas_2018",         "LUAD",
  "lusc_tcga_pan_can_atlas_2018",         "LUSC",
  "ov_tcga_pan_can_atlas_2018",           "OV",
  "blca_tcga_pan_can_atlas_2018",         "BLCA",
  "skcm_tcga_pan_can_atlas_2018",         "SKCM",
  "stad_tcga_pan_can_atlas_2018",         "STAD",
  "lihc_tcga_pan_can_atlas_2018",         "LIHC",
  "hnsc_tcga_pan_can_atlas_2018",         "HNSC",
  "lgg_tcga_pan_can_atlas_2018",          "LGG",
  "ucec_tcga_pan_can_atlas_2018",         "UCEC",
  "esca_tcga_pan_can_atlas_2018",         "ESCA",
  "gbm_tcga_pan_can_atlas_2018",          "GBM",
  "cesc_tcga_pan_can_atlas_2018",         "CESC",
  "ucs_tcga_pan_can_atlas_2018",          "UCS",
  "kirc_tcga_pan_can_atlas_2018",         "KIRC",
  "kirp_tcga_pan_can_atlas_2018",         "KIRP",
  "acc_tcga_pan_can_atlas_2018",          "ACC",
  "prad_tcga_pan_can_atlas_2018",         "PRAD",
  "thca_tcga_pan_can_atlas_2018",         "THCA",
  "tgct_tcga_pan_can_atlas_2018",         "TGCT",
  "pcpg_tcga_pan_can_atlas_2018",         "PCPG",
  "paad_tcga_pan_can_atlas_2018",         "PAAD",
  "meso_tcga_pan_can_atlas_2018",         "MESO",
  "dlbc_tcga_pan_can_atlas_2018",         "DLBC",
  "chol_tcga_pan_can_atlas_2018",         "CHOL",
  "kich_tcga_pan_can_atlas_2018",         "KICH",
  "uvm_tcga_pan_can_atlas_2018",          "UVM"
)

fusions <- fusions %>%
  left_join(study_map, by = "studyId") %>%
  mutate(
    cancer_type  = coalesce(cancer_type,
                            toupper(gsub("_tcga.*", "", studyId))),
    # Gene-pair name (directional: 5' gene – 3' gene)
    fusion_name  = paste0(site1HugoSymbol, "–", site2HugoSymbol),
    # Expression proxy: RNA split reads supporting the breakpoint
    expression   = as.numeric(tumorSplitReadCount),
    frame_status = case_when(
      site2EffectOnFrame == "in-frame"    ~ "In-frame",
      site2EffectOnFrame == "frameshift"  ~ "Frameshift",
      TRUE                                ~ NA_character_
    )
  )


clinical_clean <- clinical %>%
  rename(patientId = submitter_id) %>%
  mutate(
    os_event   = as.integer(vital_status == "Dead"),
    # Use days_to_death if dead; otherwise days_to_last_follow_up
    os_time    = as.numeric(
      if_else(vital_status == "Dead",
              as.numeric(days_to_death),
              as.numeric(days_to_last_follow_up))
    ),
    os_time_mo = os_time / 30.4375,            # days → months
    age_years  = as.numeric(age_at_diagnosis) / 365.25   # days → years
  ) %>%
  filter(!is.na(os_time), os_time > 0)



# Fusion burden per patient (and per sample)
burden_pat  <- fusions %>%
  count(patientId, cancer_type, name = "n_fusions")

burden_samp <- fusions %>%
  count(sampleId, cancer_type, name = "n_fusions")




p1a <- burden_pat %>%
  group_by(cancer_type) %>%
  summarise(med = median(n_fusions), n_pts = n(), .groups = "drop") %>%
  arrange(desc(med)) %>%
  mutate(cancer_type = fct_inorder(cancer_type)) %>%
  ggplot(aes(cancer_type, med, fill = cancer_type)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f\nn=%d", med, n_pts)),
            vjust = -0.15, size = 2.6, lineheight = 0.85) +
  scale_fill_manual(values = CANCER_COLORS, na.value = "grey70") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = " Median fusions per patient by cancer type",
       x = NULL, y = "Median fusions / patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 1b: Distribution histogram (log₁₀ scale)
med_fus <- median(burden_pat$n_fusions)
p1b <- burden_pat %>%
  ggplot(aes(n_fusions)) +
  geom_histogram(bins = 45, fill = "steelblue3", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = med_fus, colour = "red", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = med_fus * 1.8, y = Inf, vjust = 1.8,
           label = paste0("Median = ", med_fus), colour = "red", size = 3.5) +
  scale_x_log10(labels = label_comma()) +
  labs(title = "Distribution of fusions per patient (pan-cancer)",
       x = "Fusions per patient (log₁₀)", y = "Patients")

# 1c: Samples per patient (multi-sample check)
p1c <- fusions %>%
  group_by(patientId) %>%
  summarise(n_samples = n_distinct(sampleId), .groups = "drop") %>%
  count(n_samples) %>%
  ggplot(aes(factor(n_samples), n, fill = factor(n_samples))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Samples per patient", x = "Samples", y = "Patients")

p_q1 <- p1a / (p1b | p1c) +
  plot_annotation(title = "Fusion Burden Across TCGA",
                  tag_levels = "a")
save_plot(p_q1, "Q1_fusion_burden.png", w = 14, h = 12)


message("\nSurvival analysis ──")

surv_df <- burden_pat %>%
  inner_join(
    clinical_clean %>% select(patientId, os_event, os_time_mo),
    by = "patientId"
  ) %>%
  filter(!is.na(os_time_mo), !is.na(os_event), os_time_mo > 0) %>%
  mutate(
    fusion_group = factor(
      if_else(n_fusions > median(n_fusions), "High (>median)", "Low (≤median)"),
      levels = c("Low (≤median)", "High (>median)")
    )
  )

message("  Patients in KM analysis: ", nrow(surv_df),
        " | Deaths: ", sum(surv_df$os_event))

fit_pan <- survfit(Surv(os_time_mo, os_event) ~ fusion_group, data = surv_df)

p_km_pan <- ggsurvplot(
  fit_pan, data = surv_df,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  palette           = c("#1565C0", "#C62828"),
  legend.labs       = c("Low fusion burden", "High fusion burden"),
  surv.median.line  = "hv",
  xlab              = "Time (months)",
  ylab              = "Overall survival probability",
  title             = "Pan-cancer OS by fusion burden (median split)",
  ggtheme           = theme_tcga
)
save_km(p_km_pan, "Q2_KM_OS_pancancer.png")


message("\nExpression vs. frame status ──")

# Expression = log(1 + tumorSplitReadCount); filter rows with known frame
expr_df <- fusions %>%
  filter(!is.na(frame_status)) %>%
  mutate(log_expr = log1p(expression))

# Pan-cancer Wilcoxon p-value
pan_wilcox_p <- wilcox.test(log_expr ~ frame_status, data = expr_df)$p.value

p3a <- expr_df %>%
  ggplot(aes(frame_status, log_expr, fill = frame_status)) +
  geom_violin(trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, fill = "white",
               outlier.size = 0.4, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("In-frame" = "#2E7D32", "Frameshift" = "#BF360C")) +
  annotate("text", x = 1.5, y = max(expr_df$log_expr, na.rm = TRUE) * 0.97,
           label = paste0("Wilcoxon p = ", formatC(pan_wilcox_p,
                          format = "e", digits = 2)),
           size = 4.2, fontface = "bold") +
  labs(title = "Fusion expression by frame status (pan-cancer)",
       subtitle = "Expression proxy = log(1 + tumor split read count)",
       x = NULL, y = "log(1 + split reads)", fill = NULL) +
  theme(legend.position = "none")

# Per-cancer Wilcoxon p-values (for annotation)
per_cancer_wilcox <- expr_df %>%
  filter(cancer_type %in% c("BRCA", "SARC", "COADREAD", "LUAD")) %>%
  group_by(cancer_type) %>%
  filter(n_distinct(frame_status) == 2) %>%
  summarise(
    p_val   = tryCatch(wilcox.test(log_expr ~ frame_status)$p.value, error = \(e) NA_real_),
    y_pos   = max(log_expr, na.rm = TRUE) * 0.97,
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE           ~ "ns"
  ))

p3b <- expr_df %>%
  filter(cancer_type %in% c("BRCA", "SARC", "COADREAD", "LUAD")) %>%
  ggplot(aes(frame_status, log_expr, fill = frame_status)) +
  geom_violin(trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.14, fill = "white",
               outlier.size = 0.3, outlier.alpha = 0.3) +
  geom_text(data = per_cancer_wilcox,
            aes(x = 1.5, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  facet_wrap(~cancer_type, scales = "free_y") +
  scale_fill_manual(values = c("In-frame" = "#2E7D32", "Frameshift" = "#BF360C")) +
  labs(title = "Expression vs. frame status — selected cancer types",
       x = NULL, y = "log(1 + split reads)", fill = NULL) +
  theme(legend.position  = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

p_q3 <- p3a | p3b
save_plot(p_q3, "Q3_expression_vs_frame.png", w = 15, h = 6)


message("\n── Top fusions per cancer type ──")

# "Highly expressed" defined by median split reads across patients carrying that fusion
top_expressed <- fusions %>%
  group_by(cancer_type, fusion_name) %>%
  summarise(
    n_patients   = n_distinct(patientId),
    med_expr     = median(expression, na.rm = TRUE),
    pct_inframe  = mean(frame_status == "In-frame", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(cancer_type) %>%
  slice_max(med_expr, n = 10, with_ties = FALSE) %>%
  ungroup()

make_top_plot <- function(ct) {
  col <- CANCER_COLORS[ct]
  if (is.na(col)) col <- "steelblue"
  top_expressed %>%
    filter(cancer_type == ct) %>%
    arrange(med_expr) %>%
    mutate(fusion_name = fct_inorder(fusion_name)) %>%
    ggplot(aes(med_expr, fusion_name, fill = pct_inframe)) +
    geom_col() +
    geom_text(aes(label = paste0("n=", n_patients)),
              hjust = -0.15, size = 3) +
    scale_fill_gradient(low = "grey85", high = col,
                        labels = percent, name = "% in-frame",
                        limits = c(0, 1)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = ct, x = "Median split reads", y = NULL) +
    theme(legend.position = "right")
}

p4_brca <- make_top_plot("BRCA")
p4_sarc <- make_top_plot("SARC")
p4_coad <- make_top_plot("COADREAD")
p4_luad <- make_top_plot("LUAD")

p_q4 <- (p4_brca | p4_sarc) / (p4_coad | p4_luad) +
  plot_annotation(
    title    = "Top 10 highest-expressed fusions per cancer type",
    subtitle = "Bar fill = fraction of occurrences that are in-frame"
  )
save_plot(p_q4, "Q4_top_expressed_fusions.png", w = 16, h = 12)

# Also: pan-cancer recurrence heatmap (top 5 by patient count)
top_recurrent <- fusions %>%
  count(cancer_type, fusion_name, name = "n_pts") %>%
  group_by(cancer_type) %>%
  slice_max(n_pts, n = 5, with_ties = FALSE) %>%
  ungroup()

p4_heat <- top_recurrent %>%
  ggplot(aes(cancer_type, fusion_name, fill = n_pts)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n_pts), size = 2.6) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Patients") +
  labs(title = "Top 5 most recurrent fusions per cancer type",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))
save_plot(p4_heat, "Q4_recurrent_fusions_heatmap.png", w = 15, h = 10)


message("\n── Fusions vs. age ──")

age_df <- burden_pat %>%
  inner_join(clinical_clean %>% select(patientId, age_years), by = "patientId") %>%
  filter(!is.na(age_years), age_years > 0, age_years < 120)

# 5a: Scatter + lm per selected cancer
p5a <- age_df %>%
  filter(cancer_type %in% c("BRCA", "SARC", "COADREAD", "LUAD")) %>%
  ggplot(aes(age_years, n_fusions, colour = cancer_type)) +
  geom_point(alpha = 0.35, size = 1.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  facet_wrap(~cancer_type, scales = "free") +
  scale_colour_manual(values = CANCER_COLORS) +
  labs(title = "Fusions per patient vs. age (selected cancers)",
       x = "Age at diagnosis (years)", y = "Fusions per patient") +
  theme(legend.position = "none")

# 5b: Pan-cancer by age decade
p5b <- age_df %>%
  mutate(age_bin = cut(age_years,
                       breaks = c(seq(20, 80, 10), Inf),
                       labels = paste0(seq(20, 80, 10), "s"),
                       right  = FALSE)) %>%
  filter(!is.na(age_bin)) %>%
  ggplot(aes(age_bin, n_fusions)) +
  geom_boxplot(fill = "steelblue3", alpha = 0.7, outlier.size = 0.5) +
  labs(title = "Fusion burden by age decade (pan-cancer)",
       x = "Age at diagnosis", y = "Fusions per patient")

# 5c: Spearman ρ lollipop across all cancers
age_corr <- age_df %>%
  group_by(cancer_type) %>%
  filter(n() >= 10) %>%
  summarise(
    rho     = cor(age_years, n_fusions, method = "spearman", use = "complete.obs"),
    p_value = tryCatch(
      cor.test(age_years, n_fusions, method = "spearman", exact = FALSE)$p.value,
      error = \(e) NA_real_
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(rho)

message("\n  Age–fusion Spearman correlations (top hits):")
print(age_corr %>% filter(p_value < 0.05) %>% arrange(p_value), n = 20)

p5c <- age_corr %>%
  mutate(cancer_type = fct_inorder(cancer_type),
         sig = !is.na(p_value) & p_value < 0.05) %>%
  ggplot(aes(rho, cancer_type, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = rho, y = cancer_type, yend = cancer_type),
               linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                      labels = c("p ≥ 0.05", "p < 0.05"), name = NULL) +
  labs(title = "Spearman ρ (age vs. fusions) by cancer type",
       x = "Spearman ρ", y = NULL)

p_q5 <- p5a / (p5b | p5c) +
  plot_annotation(title = "Fusion Count vs. Age at Diagnosis",
                  tag_levels = "a")
save_plot(p_q5, "Q5_fusions_vs_age.png", w = 14, h = 13)


message("\n── Cancer-specific survival analyses ──")

run_cancer_km <- function(cancer_abbr, cancer_clin_df) {
  message("  Processing: ", cancer_abbr)

  f <- fusions %>% filter(cancer_type == cancer_abbr)

  clin <- cancer_clin_df %>%
    rename(patientId = submitter_id) %>%
    mutate(
      os_event   = as.integer(vital_status == "Dead"),
      os_time_mo = as.numeric(
        if_else(vital_status == "Dead",
                as.numeric(days_to_death),
                as.numeric(days_to_last_follow_up))
      ) / 30.4375
    ) %>%
    filter(!is.na(os_time_mo), os_time_mo > 0)

  b <- f %>% count(patientId, name = "n_fusions")

  sdata <- b %>%
    inner_join(clin, by = "patientId") %>%
    mutate(
      fusion_group = factor(
        if_else(n_fusions > median(n_fusions), "High (>median)", "Low (≤median)"),
        levels = c("Low (≤median)", "High (>median)")
      )
    )

  if (nrow(sdata) < 20 || n_distinct(sdata$fusion_group) < 2) {
    message("    Skipped — insufficient data"); return(invisible(NULL))
  }

  message("    n=", nrow(sdata), " | events=", sum(sdata$os_event))

  fit <- survfit(Surv(os_time_mo, os_event) ~ fusion_group, data = sdata)

  p <- ggsurvplot(
    fit, data = sdata,
    pval              = TRUE,
    pval.method       = TRUE,
    conf.int          = TRUE,
    risk.table        = TRUE,
    risk.table.height = 0.28,
    palette           = c("#1565C0", "#C62828"),
    legend.labs       = c("Low fusion burden", "High fusion burden"),
    surv.median.line  = "hv",
    xlab              = "Time (months)",
    ylab              = "OS probability",
    title             = paste0(cancer_abbr, " — OS by fusion burden",
                               "  (n=", nrow(sdata),
                               ", events=", sum(sdata$os_event), ")"),
    ggtheme           = theme_tcga
  )
  save_km(p, paste0("CS_", cancer_abbr, "_KM_OS.png"))
}

run_cancer_km("BRCA",     brca_clin)
run_cancer_km("SARC",     sarc_clin)
run_cancer_km("COADREAD", coadread_clin)
run_cancer_km("LUAD",     luad_clin)


message("\n── Frame status vs. OS ──")

frame_surv <- fusions %>%
  filter(!is.na(frame_status)) %>%
  group_by(patientId) %>%
  summarise(
    pct_inframe = mean(frame_status == "In-frame"),
    n_fusions   = n(),
    .groups = "drop"
  ) %>%
  inner_join(
    clinical_clean %>% select(patientId, os_event, os_time_mo),
    by = "patientId"
  ) %>%
  filter(!is.na(os_time_mo), os_time_mo > 0) %>%
  mutate(
    frame_group = factor(
      if_else(pct_inframe > 0.5, "Majority in-frame", "Majority frameshift"),
      levels = c("Majority frameshift", "Majority in-frame")
    )
  )

fit_frame <- survfit(Surv(os_time_mo, os_event) ~ frame_group, data = frame_surv)

p_km_frame <- ggsurvplot(
  fit_frame, data = frame_surv,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  palette           = c("#BF360C", "#2E7D32"),
  legend.labs       = c("Majority frameshift", "Majority in-frame"),
  surv.median.line  = "hv",
  xlab              = "Time (months)",
  ylab              = "OS probability",
  title             = "Pan-cancer OS by dominant fusion frame status",
  ggtheme           = theme_tcga
)
save_km(p_km_frame, "frame_status_OS.png")


message("\n── AJCC stage vs. fusion burden ──")

cancer_clin_list <- list(
  BRCA     = brca_clin,
  SARC     = sarc_clin,
  COADREAD = coadread_clin,
  LUAD     = luad_clin
)

stage_plots <- lapply(names(cancer_clin_list), function(ct) {
  clin_df <- cancer_clin_list[[ct]]

  if (!"ajcc_pathologic_stage" %in% names(clin_df)) return(NULL)

  clin <- clin_df %>%
    rename(patientId = submitter_id) %>%
    select(patientId, ajcc_pathologic_stage) %>%
    filter(!is.na(ajcc_pathologic_stage),
           !ajcc_pathologic_stage %in% c("Not Reported", "Stage X", ""))

  d <- fusions %>%
    filter(cancer_type == ct) %>%
    count(patientId, name = "n_fusions") %>%
    inner_join(clin, by = "patientId")

  if (nrow(d) < 10) return(NULL)

  d %>%
    ggplot(aes(ajcc_pathologic_stage, n_fusions,
               fill = ajcc_pathologic_stage)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) +
    scale_fill_brewer(palette = "RdYlBu") +
    labs(title = ct, x = NULL, y = "Fusions / patient") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
})

stage_plots <- Filter(Negate(is.null), stage_plots)

if (length(stage_plots) >= 2) {
  p_stage <- wrap_plots(stage_plots) +
    plot_annotation(
      title    = "Fusion burden by AJCC pathologic stage",
      subtitle = "Does more advanced disease carry more fusions?"
    )
  save_plot(p_stage, "stage_vs_fusions.png", w = 14, h = 8)
}


message("\nIn-frame fraction by cancer type ──")

inframe_summary <- fusions %>%
  filter(!is.na(frame_status)) %>%
  count(cancer_type, frame_status) %>%
  group_by(cancer_type) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

p_inframe <- inframe_summary %>%
  filter(frame_status == "In-frame") %>%
  arrange(pct) %>%
  mutate(cancer_type = fct_inorder(cancer_type)) %>%
  ggplot(aes(pct, cancer_type, fill = cancer_type)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = percent(pct, accuracy = 1)),
            hjust = -0.1, size = 3.2) +
  scale_x_continuous(labels = percent, limits = c(0, 0.85)) +
  scale_fill_manual(values = CANCER_COLORS, na.value = "grey70") +
  geom_vline(xintercept = mean(inframe_summary$pct[inframe_summary$frame_status == "In-frame"],
                                na.rm = TRUE),
             linetype = "dashed", colour = "grey30") +
  labs(title = "% in-frame fusions by cancer type",
       subtitle = "Dashed line = pan-cancer mean",
       x = "Fraction in-frame", y = NULL)
save_plot(p_inframe, "inframe_fraction_by_cancer.png", w = 9, h = 7)


message("\n", strrep("─", 55))
message("All plots saved to: ", normalizePath(OUTDIR))
message("\nFiles created:")
list.files(OUTDIR, full.names = FALSE) %>%
  sort() %>%
  paste0("  ", .) %>%
  paste(collapse = "\n") %>%
  message()
