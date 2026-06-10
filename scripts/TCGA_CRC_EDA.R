library(dplyr)
library(ggplot2)
library(tidyverse)

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
  ) %>%
  distinct(submitter_id, .keep_all = TRUE) %>%  # one row per patient
  filter(!is.na(age_at_index))                  # remove missing ages

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


# Add total fusion per age and mean fusion per age

age_table <- age_table %>%
  group_by(age_band) %>%
  mutate(
    nPatient_per_age_group = n_distinct(submitter_id),
    total_fusions = sum(nfusion_per_patient),
    mean_fusions_per_patient = total_fusions / nPatient_per_age_group
  ) %>%
  ungroup()
age_bar_plot<-ggplot(age_table,
                     aes(x=age_band,
                         y=mean_fusions_per_patient))+geom_col()+
  theme_classic()+ labs(
    x = "Age Group",
    y = "mean Fusion per Age Group"
  )



# Add scatter plot to see if there is a patttern
age_scater_plot<-ggplot(age_table,
                     aes(x=age_at_index,
                         y=nfusion_per_patient))+geom_point()+
  geom_smooth(method = lm,color="red", se=T)+
  theme_classic()+ labs(
    x = "Age",
    y = "Sum of Fusion per Patient"
  )

age_scater_plot

ggsave("age_bar_plot.png",age_scater_plot, dpi=300,
       width = 6,
       height = 4)


#==========================
# Prognosis and Fusion per age box plot
#==========================
print("Do poor responders have higher number of fusions ?")

response_table <- master %>%
  distinct(submitter_id, .keep_all = TRUE) %>%
  filter(!is.na(follow_ups_disease_response)) %>%
  filter(! follow_ups_disease_response=="Unknown")%>%
  select(
    submitter_id,
    follow_ups_disease_response,
    nfusion_per_patient
  )


unique(response_table$follow_ups_disease_response)
# Boxplot + individual patient points
ggplot(response_table,
       aes(x = follow_ups_disease_response,
           y = nfusion_per_patient),
       box.colour = c("pink","purple")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  theme_classic() +
  labs(
    x = "Disease response group",
    y = "Number of fusions per patient",
    title = "Fusion burden by disease response"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



print("Is higher number of fusion correlated to poor response")




# Do higher number of fusions lead to poorer prognosis?----