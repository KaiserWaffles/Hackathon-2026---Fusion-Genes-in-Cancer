import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
fusion_df = pd.read_csv(r"data/cbioportalR_TCGA_fusions.csv")
#print(fusion_df.columns)

clinical = pd.read_csv(r"data/clinical_shared_columns.csv")
print(fusion_df ['patientId'].head())
print(clinical ['submitter_id'].head())

#fusion burden per patient
fusion_df['fusion_name'] = (
    fusion_df['site1HugoSymbol']
    + '--' +
    fusion_df['site2HugoSymbol']
)
fusion_per_patient = (
    fusion_df[['patientId', 'fusion_name']]
    .drop_duplicates()
    .groupby('patientId')
    .size()
    .reset_index(name='n_fusions')
)
merged = fusion_per_patient.merge(
    clinical,
    left_on='patientId',
    right_on='submitter_id',
    how='inner'
)

merged['time'] = merged['days_to_death'].fillna(merged['days_to_last_follow_up'])
merged['event'] = (merged['vital_status'] == 'Dead').astype(int)
merged = merged.dropna(subset=['time', 'n_fusions'])
merged = merged.dropna(subset=['ajcc_pathologic_stage'])

plt.figure(figsize=(10,6))

sns.boxplot(
    data=merged,
    x='ajcc_pathologic_stage',
    y='n_fusions'
)

plt.xticks(rotation=45)
plt.title("Fusion Burden across Cancer Stages")

plt.savefig("stage_vs_fusions.png", dpi=300)
plt.show()