library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# Set working directory
setwd("D:/gene_fusion")


# Load the data----
crc.metadata <- read.csv("D:/gene_fusion/data/TCGA/TCGA/clinical/TCGA-COAD_clinical.csv")
crc.fusion.info <- read.csv("D:/gene_fusion/data/TCGA/TCGA/cbioportalR_TCGA_fusions.csv")

# Get number of fusions per sample & per patient----
# Per sample
crc.fusion.info <- crc.fusion.info %>%
  group_by(sampleId) %>%
  mutate(nfusion_per_sample = n()) %>%
  ungroup()

# per patient
crc.fusion.info <- crc.fusion.info%>%
  group_by(patientId)%>%
  mutate(nfusion_per_patient =n())%>%
  ungroup()

# how many sample per patient
crc.fusion.info <- crc.fusion.info %>%
  group_by(patientId) %>%
  mutate(nsample_per_patient = n_distinct(sampleId)) %>%
  ungroup()

unique(crc.fusion.info$nsample_per_patient) # all patient have given one sample s it is the same


#========================================
# Colorectal Cancer EDA===================
#=============================================

# Extract only colorectal cancer patients fusion events----
crc.fusion.info <- crc.fusion.info[crc.fusion.info$patientId 
                                   %in% crc.metadata$submitter_id,]



# Merge two together 
master <- merge(
  x = crc.metadata,
  y = crc.fusion.info,
  by.x = "submitter_id",
  by.y = "patientId"
)


#=========================
# Fusion number per patient Density plot
#=========================

fusion.density.p<-ggplot(crc.fusion.info,
       aes(x = nfusion_per_patient)) +
  geom_density()+
  theme_classic() +
  labs(
    x = "Number of Fusions",
    y = "Probability Density",
    title = "Fusion Density"
  )


ggsave("fusion.density.p.png",fusion.density.p, dpi=300,
       width = 6,
       height = 4)

#===========================================
# Age=====
#===========================================

# I will creat bands of age ----

# add age bands
age_table <- master %>%
  select(
    submitter_id,
    age_at_index,
    age_at_diagnosis,
    nfusion_per_patient
  ) %>% distinct(submitter_id, .keep_all = TRUE) %>%  
  filter(!is.na(age_at_index))                

# Create 10-year age bands
age_breaks <- seq(
  from = floor(min(age_table$age_at_index)),
  to   = ceiling(max(age_table$age_at_index)) + 10,
  by   = 10
)

age_table <- age_table %>%
  mutate(
    age_band = cut(
      age_at_index,
      breaks = age_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  )



#= Boxplot====
ggplot(age_table,
       aes(age_band,
           nfusion_per_patient))+
  geom_violin(trim=FALSE)+
  geom_boxplot()+
  theme_classic()

#==========================
# Prognosis and Fusion per age box plot
#==========================
print("Do poor responders have higher number of fusions ?")

response_table <- master %>%
  distinct(submitter_id, .keep_all = TRUE) %>%
  filter(!is.na(follow_ups_disease_response)) %>%
  filter(! follow_ups_disease_response=="Unknown")%>%
  dplyr:: select(
    submitter_id,
    follow_ups_disease_response,
    nfusion_per_patient
  )
response_table$clean_response_cat <- ifelse(
response_table$follow_ups_disease_response== "TF-Tumor Free","Responders","Non-Responders")

# Boxplot + individual patient points-----
p<-ggplot(response_table,
       aes(x = clean_response_cat,
           y = nfusion_per_patient)) +
  geom_boxplot(outlier.shape = NA,,color = c("red","blue")) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +

  theme_classic(base_size = 13, base_family = "Arial") +
  labs(
    x = NULL,
    y = "Number of fusions per patient",
    title = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0),
    legend.position = "none"
  )


my_comparisons<- c("Responders","Non-Responders")
p.response<-p+stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",
  label = "p.format"
)+  stat_compare_means(label.y = 13)    

ggsave("plot/response.boxplot.png",
       p.response,
       dpi=300,
       width = 6,
       height = 5)



#stats----
response_table %>%
  group_by(clean_response_cat) %>%
  summarise(
    n=n(),
    mean=mean(nfusion_per_patient),
    median=median(nfusion_per_patient),
    sd=sd(nfusion_per_patient),
    min=min(nfusion_per_patient),
    max=max(nfusion_per_patient)
  )
