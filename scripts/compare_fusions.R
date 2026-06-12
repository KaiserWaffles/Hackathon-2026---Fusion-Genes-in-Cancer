#rm(list = ls())


library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggvenn)        
library(stringr)
library(gprofiler2)

setwd("D:/gene_fusion")
#########
# Load data####
######
short.og <- read.csv("all_fusions_collated.csv")



long<- readxl::read_xlsx("data/CRC_human.xlsx",
                               sheet = "Table S6")
names(long)<-long[1,]
long<-long[-1,]
colnames(long)
# only keep gene
long <- long[,c("Gene1","Gene2")]
long <- long%>%rename(
    gene1            = Gene1,
    gene2            = Gene2
  )

dim(long)
############################################
# Uppercase for mouse gene
#==========================================
#gene 1
mouse_gene1 <- short.og[short.og$dataset=="mouse_ifn",]%>% pull (gene1_clean)
ortholog_gene_1 <- gorth(
  query = mouse_gene1,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)%>% select(c(input,ortholog_name))

colnames(ortholog_gene_1)<- c("mouse_gene1","h_gene1")

#gene2
mouse_gene2 <- short.og[short.og$dataset=="mouse_ifn",]%>%
  pull (gene2_clean)

ortholog_gene_2 <- gorth(
  query = mouse_gene2,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)%>% select(c(input,ortholog_name))
colnames(ortholog_gene_2)<- c("mouse_gene2","h_gene2")



# replace gene1 with human ortholog
mouse_db <- short.og[short.og$dataset=="mouse_ifn",]%>%
  select(gene1_clean,gene2_clean,confidence,dataset)
dim(mouse_db)
dim(short.og)
mouse_db_1 <- merge(
  x = mouse_db,
  y = ortholog_gene_1,
  by.x = "gene1_clean",
  by.y = "mouse_gene1"
)

# replace gene2 
mouse_db_human <- merge(
  x = mouse_db_1,
  y = ortholog_gene_2,
  by.x = "gene2_clean",
  by.y = "mouse_gene2"
)

# keep only human gene names
colnames(mouse_db_human)
mouse_db_human <- mouse_db_human %>%
  select(
    h_gene1, h_gene2,confidence,dataset
  )

head(mouse_db_human)
colnames(mouse_db_human)<-c("gene1","gene2","confidence","dataset"   )


#================
# Create short clean
#==============
short<-short.og[!short.og$dataset=="mouse_ifn",]%>%select(
  gene1_clean,gene2_clean,confidence,dataset
)
unique(short$dataset)
dim(short) #458

colnames(short)<-c("gene1","gene2","confidence","dataset")
short.clean <-rbind(short,mouse_db_human) 
unique(short.clean$dataset)
dim(short.clean) #35138     

#Alphabetically sort fusion keys to collapse direction
short.clean <- short.clean %>%
  mutate(fusion_key = paste(pmin(gene1, gene2),
                            pmax(gene1, gene2),
                            sep = "_"))
length(unique(short.clean$fusion_key))#339
long <- long %>%
  mutate(fusion_key = paste(pmin(gene1, gene2),
                            pmax(gene1, gene2),
                            sep = "_"))
length(unique(long$fusion_key))#1116

# Keep only high confidance short reads 
"short_filtered <- short %>%
  filter(
    confidence    == "high"
  )"

# How many fusions remain per study after filtering?
"short_filtered %>%
  count(condition, name = "n_high_confidence") %>%
  print()
"
#==================================
# Each study alone ====
#======================
mouse_ifn  <- short.clean %>%
  filter(dataset=="mouse_ifn") %>%
  pull(fusion_key) %>% unique()

human_bile   <- short.clean %>%
  filter(dataset== "human_bile") %>%
  pull(fusion_key) %>% unique()

human_antiEGFR   <- short.clean %>%
  filter(dataset== "human_antiEGFR") %>%
  pull(fusion_key) %>% unique()



#===================
# Ven####
#=====================
keys_short_all   <- unique(short.clean$fusion_key)
keys_long        <- unique(long$fusion_key)

shared_all     <- intersect(keys_short_all, keys_long)
short_only     <- setdiff(keys_short_all, keys_long)
long_only      <- setdiff(keys_long, keys_short_all)

"keys_short_all_gene <-unique(short.clean$gene1)
keys_long_gene        <- unique(long$gene1)
shared_all_gene    <- intersect(keys_short_all_gene,
                                keys_long_gene)


keys_short_all_gene2 <-unique(short.clean$gene2)
keys_long_gene2        <- unique(long$gene2)
shared_all_gene2    <- intersect(keys_short_all_gene2,
                                keys_long_gene2)
"

# Venn diagram — short reads vs long reads
venn_data <- list(
  `Short reads` = keys_short_all,
  `Long reads` = keys_long
)

ggvenn(venn_data, fill_color = c("#2E75B6", "#70AD47"),
       stroke_size = 0.5, set_name_size = 4) +
  labs(title = "Fusion overlap:") +
  theme_void(base_size = 12)

ggsave("plot/venn_platform_overlap.pdf",dpi=300, width = 6, height = 5)

# 4c. Three-way Venn across studies (human organoids only — same species as long reads)
venn_human <- list(
  `Mouse IFN` = mouse_ifn,
  `Human bile` = human_bile,
  `Human antiEGFR`    = human_antiEGFR
)





ggvenn(venn_human, fill_color = c("#2E75B6", "#ED7D31", "#70AD47"),
       stroke_size = 0.5, set_name_size = 3.5) +
  labs(title = "Fusion overlap:") +
  theme_void(base_size = 12)

ggsave("venn_human_three_way.pdf", width = 7, height = 6)


shared_human_bile_antiegfr   <- intersect(human_bile,
                                 human_antiEGFR)
