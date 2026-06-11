library(tidyverse)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db )
library(enrichplot)

# Set working directory
setwd("D:/gene_fusion")


#================
# Read samples 
#================
ss<- readxl::read_xlsx("data/bulk_RNA/mouse/ss.xlsx")

data_path= "data/bulk_RNA/mouse"
# read files
sample_names <- ss$sample
samples_list <- list()
for (n in sample_names){
  
  sample <- read_tsv(file.path(data_path,
                               paste0(n,".tsv")))
  
  samples_list[[n]]<- sample
}


# check if all sample have ssample 

all_same_cols <- all(
  sapply(samples_list, function(x)
    identical(colnames(x), colnames(samples_list[[1]])))
)

all_same_cols # T




# Safe merge 
merged_samples <- bind_rows(samples_list)




#===================
# go
#==============
#control----
fusion_genes <- c(merged_samples$`#gene1`) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enriched_gene1 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Mm.eg.db, ont = "MF",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

p_gene1 <- enrichplot::dotplot(enriched_gene1)
p_gene1


# control gene 2----
fusion_genes <- c(merged_samples$gene2) %>% unique()
markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enriched_gene2 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Mm.eg.db, ont = "MF",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

p_gene2 <- enrichplot::dotplot(enriched_gene2)
p_gene2
ggsave("plot/CRC_mouse_GO_gene2.png",
       p_gene2,
       dpi=300,
       width = 8,
       height = 6)


# ALL genes
fusion_genes <- c(merged_samples$gene2,merged_samples$`#gene1`) %>% unique()
markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
enriched = enrichGO(gene = markers$ENTREZID,OrgDb = org.Mm.eg.db, ont = "MF",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

p <- enrichplot::dotplot(enriched)
p
ggsave("plot/CRC_mouse_GO.png",
       p,
       dpi=300,
       width = 8,
       height = 6)

