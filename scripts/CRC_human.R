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
# Read xslx
#=======================
crc_fusion<- readxl::read_xlsx("data/CRC_human.xlsx",
                               sheet = "Table S6")
names(crc_fusion)<-crc_fusion[1,]
crc_fusion<-crc_fusion[-1,]


# =====================
# What gene molecular functions do these fusion genes have
#=======================
# gene 1----

fusion_genes <- c(crc_fusion$Gene1) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched_ctrl_gene1 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE)

p_crl_gene1 <- enrichplot::dotplot(enriched_ctrl_gene1)
p_crl_gene1
ggsave("plot/CRC_human_gene1_GO.png",
       p_crl_gene1,
       dpi=300,
       width = 8,
       height = 6)


# genen 2----


fusion_genes <- c(crc_fusion$Gene2) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched_ctrl_gene2 = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE)

p_crl_gene2 <- enrichplot::dotplot(enriched_ctrl_gene2)
p_crl_gene2
ggsave("plot/CRC_human_gene2_GO.png",
       p_crl_gene2,
       dpi=300,
       width = 8,
       height = 6)

# gene 1 and gene 2-----

fusion_genes <- c(crc_fusion$Gene2,
                  crc_fusion$Gene2) %>% unique()

markers = bitr(fusion_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enriched = enrichGO(gene = markers$ENTREZID,OrgDb = org.Hs.eg.db, ont = "MF",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE)

p <- enrichplot::dotplot(enriched)
p
ggsave("plot/CRC_human_GO.png",
       p,
       dpi=300,
       width = 8,
       height = 6)

#=============================
# mouse data ====
#===========================

