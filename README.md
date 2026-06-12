# Hackathon-2026---Fusion-Genes-in-Cancer
 In this project, we aim to explore RNA fusion events in different cancer RNA-Seq data generated at the CRUK Scotland institute

## Exploring the TCGA Data Set for CRC

Exploring the relationship between fusion burden, age, and response in Colon Adenocarcinoma dataset within The Cancer Genome Atlas (TCGA

**TCGA_CRC_EDA.R:** [TCGA_CRC_EDA.R](https://github.com/KaiserWaffles/Hackathon-2026---Fusion-Genes-in-Cancer/blob/main/scripts/TCGA_CRC_EDA.R) -Anya

## GO Ontology 
Exploring the molacular function of the genes involved in fusion using Sun, Qiang et al.(2023) for human CRC, and  Jing, Zhi-Liang et al. Scientific reports (2024) mouse data.
**Code:** [CRC_human:](https://github.com/KaiserWaffles/Hackathon-2026---Fusion-Genes-in-Cancer/blob/main/scripts/CRC_human.R) -Anya
**Code**[CRC_mouse:](https://github.com/KaiserWaffles/Hackathon-2026---Fusion-Genes-in-Cancer/blob/main/scripts/CRC_mouse.R)- Anya

## Shared Fusions Between Long Reads and Short Reads
Investigating common fusions between Long reads obtained from  [Sun, Qiang et al.(2023)] (doi:10.1186/s13073-023-01226-y) and three paired end bulk RNA seq reads that were processed through Arriba to find the fusion events. These three datasets are:
| Dataset                              | Citation                                     |
| ------------------------------------ | -------------------------------------------- |
| PRJNA1269539_arriba_human_bile_salts | **Lässle et al., 2025**                      |
| PRJNA1062304_arriba_mouse_ifn        | **Jing et al., 2024**                        |
| PRJNA1462607_arriba_human_antiEGFR   | **Unpublished dataset (PRJNA1462607, 2026)** |
**Code** [Compare_fusions.R](https://github.com/KaiserWaffles/Hackathon-2026---Fusion-Genes-in-Cancer/blob/main/scripts/compare_fusions.R)- Anya


