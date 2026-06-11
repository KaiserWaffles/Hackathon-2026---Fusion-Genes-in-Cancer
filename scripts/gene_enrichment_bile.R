library(tidyverse)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(org.Mm.eg.db )
library(enrichplot)

# Set working directory
setwd("D:/gene_fusion")

#======================
# Read TSV files for bile salt bulk rna----
#=======================
ss<- readxl::read_xlsx("data/bulk_RNA/bile_salt/ss.xlsx")

data_path= "data/bulk_RNA/bile_salt"
# read files
sample_names <- ss$sample
samples_list <- list()
for (n in sample_names){
  
  sample <- read_tsv(file.path(data_path,
                               paste0(n,".tsv")))
  
  sample$group <- ss$group[ss$sample == n]
  samples_list[[n]]<- sample
}


# check if all sample have ssample 

all_same_cols <- all(
  sapply(samples_list, function(x)
    identical(colnames(x), colnames(samples_list[[1]])))
)

all_same_cols # T


# delete breakpoint.1, and 2 becuase some tables are all NA

samples_list <- lapply(samples_list, function(x) {
  x %>%
    select(-breakpoint1, -breakpoint2)
})

# Safe merge 
merged_samples <- bind_rows(samples_list)


# GO on control----

#===================
# CONTROL-----
#==============
# Control gene 1
ctrl <- merged_samples[merged_samples$group=="control",]
fusion_genes <- c(ctrl$`#gene1`) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched_ctrl_gene1 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                    pvalueCutoff  = 0.25,
                    qvalueCutoff  = 0.25,
                    readable      = TRUE)

p_crl_gene1 <- enrichplot::dotplot(enriched_ctrl_gene1)
p_crl_gene1


# control gene 2
fusion_genes <- c(ctrl$gene2) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched_ctrl_gene2 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                               pvalueCutoff  = 0.25,
                               qvalueCutoff  = 0.25,
                               readable      = TRUE)

p_crl_gene2 <- enrichplot::dotplot(enriched_ctrl_gene2)
p_crl_gene2



#===================
# every bile acid-----
#==============
# Control gene 1
intervention <- merged_samples[merged_samples$group %in% c("CA",
                                                   "GCA",
                                                   "TCA"),]

fusion_genes <- c(intervention$`#gene1`) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched_gene1 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                               pvalueCutoff  = 0.25,
                               qvalueCutoff  = 0.25,
                               readable      = TRUE)

p_gene1 <- enrichplot::dotplot(enriched_gene1)
p_gene1


# control gene 2
fusion_genes <- c(intervention$gene2) %>% unique()
??enrichGO
markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched_gene2 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)

p_gene2 <- enrichplot::dotplot(enriched_gene2)
p_gene2

