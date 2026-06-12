#Comparing the fusion genes found between different groups
# Load required libraries
library(dplyr)
library(tidyr)
library(UpSetR)
library(stringr)
library(purrr)
library(VennDiagram)
library(gridExtra)
library(ggVennDiagram)


#Create function to separate >1 genes
gene_separator <- function(genes) {
  gene_list <- str_extract_all(genes, "[A-Za-z0-9]+(?=\\()")
  result <- unlist(gene_list)
  if(length(result) == 0) {
    result <- trimws(unlist(strsplit(genes, ",")))
  }
  return(result)
}

#Assign original TSVs to variable
mouse_files <- list.files("results/mouse_ifn_results", pattern = "\\.tsv$", full.names = TRUE)
human_bile_files <- list.files("results/human_bile_results", pattern = "\\.tsv$", full.names = TRUE)
human_egfr_files <- list.files("results/human_antiEGFR_results", pattern = "\\.tsv$", full.names = TRUE)

# Create function to extract fusions from a single file
extract_fusions <- function(file_path, sample_name) {
  df <- read.delim(file_path, sep = '\t', check.names = FALSE)
  
  # Get the fusion_genes column
  fusion_col <- grep("fusion_genes", names(df), value = TRUE)[1]
  
  fusions <- gene_separator(paste(df[[fusion_col]][df[[fusion_col]] != ""], collapse = ","))
  
  # Return as data frame
  data.frame(
    sample = sample_name,
    fusion = fusions,
    stringsAsFactors = FALSE
  )
}

# Process anti-EGFR files
antiEGFR_data <- bind_rows(lapply(human_egfr_files, function(file) {
  if(grepl("ctrl", file, ignore.case = TRUE)) {
    # Extract number from ctrl_1, ctrl_2, etc.
    rep_num <- str_extract(file, "ctrl_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("EGFR_Ctrl_", rep_num)
  } else if(grepl("antiEGFR", file, ignore.case = TRUE)) {
    # Extract number from antiEGFR_1, antiEGFR_2, etc.
    rep_num <- str_extract(file, "antiEGFR_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("EGFR_Treat_", rep_num)
  } else {
    # Fallback for any other pattern
    sample_name <- paste0("EGFR_", basename(file))
  }
  extract_fusions(file, sample_name)
}))

#Process bile data
bile_data <- bind_rows(lapply(human_bile_files, function(file) {
  basename_file <- basename(file)
  
  if(grepl("ctrl", basename_file, ignore.case = TRUE)) {
    rep_num <- str_extract(basename_file, "ctrl_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("Bile_Ctrl_", rep_num)
  } else if(grepl("_CA_", basename_file)) {
    rep_num <- str_extract(basename_file, "CA_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("Bile_CA_", rep_num)
  } else if(grepl("_GCA_", basename_file)) {
    rep_num <- str_extract(basename_file, "GCA_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("Bile_GCA_", rep_num)
  } else if(grepl("_TCA_", basename_file)) {
    rep_num <- str_extract(basename_file, "TCA_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("Bile_TCA_", rep_num)
  } else {
    sample_name <- paste0("Bile_", basename_file)
  }
  extract_fusions(file, sample_name)
}))

# Process mouse files
mouse_data <- bind_rows(lapply(mouse_files, function(file) {
  if(grepl("ctrl", file)) {
    rep_num <- str_extract(file, "ctrl_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("Mouse_Ctrl_", rep_num)
  } else {
    rep_num <- str_extract(file, "IFN_([0-9])") %>% str_extract("[0-9]")
    sample_name <- paste0("Mouse_IFN_", rep_num)
  }
  extract_fusions(file, sample_name)
}))

all_fusions <- bind_rows(antiEGFR_data, bile_data, mouse_data)

#Print out summary stats
print("Samples in dataset:")
print(unique(all_fusions$sample))
print(paste("Total fusions:", nrow(all_fusions)))
print(paste("Unique fusions:", n_distinct(all_fusions$fusion)))
print("Fusions per sample:")
print(table(all_fusions$sample))

#===============================================================================
# Group by treatment type
grouped_fusions <- all_fusions %>%
  mutate(group = case_when(
    grepl("EGFR_Ctrl", sample) ~ "EGFR_Control",
    grepl("EGFR_Treat", sample) ~ "EGFR_Treatment",
    grepl("Bile_Ctrl", sample) ~ "Bile_Control",
    grepl("Bile_CA", sample) ~ "Bile_CA",
    grepl("Bile_GCA", sample) ~ "Bile_GCA",
    grepl("Bile_TCA", sample) ~ "Bile_TCA",
    grepl("Mouse_Ctrl", sample) ~ "Mouse_Control",
    grepl("Mouse_IFN", sample) ~ "Mouse_IFN"
  )) %>%
  group_by(fusion, group) %>%
  summarise(present = 1, .groups = "drop") %>%
  distinct()

# Create grouped matrix
group_matrix <- grouped_fusions %>%
  pivot_wider(id_cols = fusion, names_from = group, values_from = present, values_fill = 0)

#==================================
# Checkpoints
treatment_fusions <- all_fusions %>%
  filter(grepl("EGFR_Treat|Bile_CA|Bile_GCA|Bile_TCA|Mouse_IFN", sample)) %>%
  mutate(group = case_when(
    grepl("EGFR_Treat", sample) ~ "EGFR_Treatment",
    grepl("Bile_CA", sample) ~ "Bile_CA",
    grepl("Bile_GCA", sample) ~ "Bile_GCA",
    grepl("Bile_TCA", sample) ~ "Bile_TCA",
    grepl("Mouse_IFN", sample) ~ "Mouse_IFN"
  ))

# Check how many groups actually have data
print(table(treatment_fusions$group))

# Create matrix and check dimensions
treatment_matrix <- treatment_fusions %>%
  distinct(group, fusion) %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = fusion, names_from = group, values_from = present, values_fill = 0)

print(paste("Matrix dimensions:", dim(treatment_matrix)))
print(paste("Number of groups (excluding fusion column):", ncol(treatment_matrix) - 1))

# Calculate how many fusions are shared between treatments
shared_counts <- treatment_fusions %>%
  group_by(fusion) %>%
  summarise(n_treatments = n_distinct(group)) %>%
  group_by(n_treatments) %>%
  summarise(count = n())

# Bar plot
ggplot(shared_counts, aes(x = factor(n_treatments), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Fusions Shared Across Treatments",
       x = "Number of Treatments Sharing the Fusion",
       y = "Number of Fusions") +
  theme_minimal()

# ============================================
# HUMAN VS MOUSE - ALL SAMPLES

# Get all human fusions (EGFR + Bile experiments)
human_fusions <- all_fusions %>%
  filter(grepl("EGFR|Bile", sample)) %>%
  pull(fusion) %>%
  unique()

# Get all mouse fusions
mouse_fusions <- all_fusions %>%
  filter(grepl("Mouse", sample)) %>%
  pull(fusion) %>%
  unique()

# Create Venn diagram
venn_human_mouse <- venn.diagram(
  x = list(
    Human = human_fusions,
    Mouse = mouse_fusions
  ),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "coral"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05),
  main = "Human vs Mouse Fusions (All Samples)",
  main.cex = 1.3
)

grid.draw(venn_human_mouse)

# Print statistics
cat("\n=== Human vs Mouse Comparison ===\n")
cat("Human-specific fusions:", length(setdiff(human_fusions, mouse_fusions)), "\n")
cat("Mouse-specific fusions:", length(setdiff(mouse_fusions, human_fusions)), "\n")
cat("Shared fusions (conserved):", length(intersect(human_fusions, mouse_fusions)), "\n")

#NO COMMON FUSION GENES BETWEEN HUMAN AND MOUSE DATA!!!

# ============================================
# ANTI-EGFR VS BILE ACIDS - ALL SAMPLES

# Get all Anti-EGFR fusions (control + treatment)
antiEGFR_all <- all_fusions %>%
  filter(grepl("EGFR", sample)) %>%
  pull(fusion) %>%
  unique()

# Get all Bile Acids fusions (control + treatment)
bile_all <- all_fusions %>%
  filter(grepl("Bile", sample)) %>%
  pull(fusion) %>%
  unique()

# Create Venn diagram
venn_human_groups <- venn.diagram(
  x = list(
    Anti_EGFR = antiEGFR_all,
    Bile_Acids = bile_all
  ),
  filename = NULL,
  col = "transparent",
  fill = c("steelblue", "forestgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05),
  main = "Anti-EGFR vs Bile Acids (All Samples)",
  main.cex = 1.3
)

grid.draw(venn_human_groups)

# Print statistics
cat("\n=== Human Group Comparison ===\n")
cat("Anti-EGFR specific fusions:", length(setdiff(antiEGFR_all, bile_all)), "\n")
cat("Bile Acids specific fusions:", length(setdiff(bile_all, antiEGFR_all)), "\n")
cat("Shared fusions:", length(intersect(antiEGFR_all, bile_all)), "\n")

shared_antiEGFR_bile <- intersect(antiEGFR_all, bile_all)

# View the shared fusions
print("Shared fusions between Anti-EGFR and Bile:")
print(shared_antiEGFR_bile)

#ONLY 1 SHARED GENE BETWEEN THE TWO GROUPS <- WASH6P :: WASH6P


# ============================================
# ANTI-EGFR: TREATMENT VS CONTROL

# Get Anti-EGFR fusions
antiEGFR_ctrl <- all_fusions %>%
  filter(grepl("EGFR_Ctrl", sample)) %>%
  pull(fusion) %>%
  unique()

antiEGFR_treat <- all_fusions %>%
  filter(grepl("EGFR_Treat", sample)) %>%
  pull(fusion) %>%
  unique()

# Create Venn diagram
venn_antiEGFR <- venn.diagram(
  x = list(
    Control = antiEGFR_ctrl,
    Treatment = antiEGFR_treat
  ),
  filename = NULL,
  col = "transparent",
  fill = c("lightblue", "steelblue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.default.pos = "outer",
  main = "Anti-EGFR: Control vs Treatment",
  main.cex = 1.2
)

grid.draw(venn_antiEGFR)

# Statistics
cat("\n=== Anti-EGFR ===\n")
cat("Control-specific:", length(setdiff(antiEGFR_ctrl, antiEGFR_treat)), "\n")
cat("Treatment-specific:", length(setdiff(antiEGFR_treat, antiEGFR_ctrl)), "\n")
cat("Shared:", length(intersect(antiEGFR_ctrl, antiEGFR_treat)), "\n")

# > cat("Control-specific:", length(setdiff(antiEGFR_ctrl, antiEGFR_treat)), "\n")
# Control-specific: 37 
# > cat("Treatment-specific:", length(setdiff(antiEGFR_treat, antiEGFR_ctrl)), "\n")
# Treatment-specific: 21 
# > cat("Shared:", length(intersect(antiEGFR_ctrl, antiEGFR_treat)), "\n")
# Shared: 20 


# ============================================
# BILE ACIDS: CONTROL VS ALL TREATMENTS

# Get Bile fusions
bile_ctrl <- all_fusions %>%
  filter(grepl("Bile_Ctrl", sample)) %>%
  pull(fusion) %>%
  unique()

bile_treat <- all_fusions %>%
  filter(grepl("Bile_CA|Bile_GCA|Bile_TCA", sample)) %>%
  pull(fusion) %>%
  unique()

# Create Venn diagram
venn_bile <- venn.diagram(
  x = list(
    Control = bile_ctrl,
    Treatment = bile_treat
  ),
  filename = NULL,
  col = "transparent",
  fill = c("lightgreen", "darkgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Bile Acids: Control vs All Treatments",
  main.cex = 1.2
)

grid.draw(venn_bile)

cat("\n=== Bile Acids ===\n")
cat("Control-specific:", length(setdiff(bile_ctrl, bile_treat)), "\n")
cat("Treatment-specific:", length(setdiff(bile_treat, bile_ctrl)), "\n")
cat("Shared:", length(intersect(bile_ctrl, bile_treat)), "\n")


# > cat("Control-specific:", length(setdiff(bile_ctrl, bile_treat)), "\n")
# Control-specific: 25 
# > cat("Treatment-specific:", length(setdiff(bile_treat, bile_ctrl)), "\n")
# Treatment-specific: 74 
# > cat("Shared:", length(intersect(bile_ctrl, bile_treat)), "\n")
# Shared: 20 


# ============================================
# MOUSE IFN: TREATMENT VS CONTROL

# Get Mouse fusions
mouse_ctrl <- all_fusions %>%
  filter(grepl("Mouse_Ctrl", sample)) %>%
  pull(fusion) %>%
  unique()

mouse_treat <- all_fusions %>%
  filter(grepl("Mouse_IFN", sample)) %>%
  pull(fusion) %>%
  unique()

# Create Venn diagram
venn_mouse <- venn.diagram(
  x = list(
    Control = mouse_ctrl,
    Treatment = mouse_treat
  ),
  filename = NULL,
  col = "transparent",
  fill = c("lightcoral", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Mouse IFN: Control vs Treatment",
  main.cex = 1.2
)

grid.draw(venn_mouse)

cat("\n=== Mouse IFN ===\n")
cat("Control-specific:", length(setdiff(mouse_ctrl, mouse_treat)), "\n")
cat("Treatment-specific:", length(setdiff(mouse_treat, mouse_ctrl)), "\n")
cat("Shared:", length(intersect(mouse_ctrl, mouse_treat)), "\n")

# > cat("Control-specific:", length(setdiff(mouse_ctrl, mouse_treat)), "\n")
# Control-specific: 96 
# > cat("Treatment-specific:", length(setdiff(mouse_treat, mouse_ctrl)), "\n")
# Treatment-specific: 83 
# > cat("Shared:", length(intersect(mouse_ctrl, mouse_treat)), "\n")
# Shared: 118 

#================================================================================
#NICER PLOTS USING GGVENNDIAGRAM

p1 = ggVennDiagram(list(Control = antiEGFR_ctrl, Treatment = antiEGFR_treat)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Anti-EGFR") +
  theme(plot.title = element_text(hjust = 0.5))

p2 = ggVennDiagram(list(Control = bile_ctrl, Treatment = bile_treat)) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(title = "Bile Acids") +
  theme(plot.title = element_text(hjust = 0.5))

p3 = ggVennDiagram(list(Control = mouse_ctrl, Treatment = mouse_treat)) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Mouse IFN") +
  theme(plot.title = element_text(hjust = 0.5))

p1
p2
p3