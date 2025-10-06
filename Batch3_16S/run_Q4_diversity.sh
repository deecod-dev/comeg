#!/bin/bash
# Q-4: Alpha diversity and evenness analysis for Batch3

# Input files
ABUNDANCE=species_rel_abundance.tsv
METADATA=Batch3.csv

# Run Python analysis
python3 << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
rel_ab = pd.read_csv("species_rel_abundance.tsv", sep="\t", index_col=0)
meta = pd.read_csv("Batch3.csv")

# Ensure sample names align
rel_ab = rel_ab.loc[:, rel_ab.columns.isin(meta['sample_id'])]

# Compute Shannon and Pielou indices
diversity = []
for sample in rel_ab.columns:
    p = rel_ab[sample].values
    p = p[p > 0]  # remove zeros
    H = -np.sum(p * np.log(p))
    S = len(p)
    J = H / np.log(S) if S > 1 else 0
    diversity.append({"sample_id": sample, "Shannon": H, "Pielou": J})

div = pd.DataFrame(diversity)
meta_div = pd.merge(meta, div, on="sample_id")

# Save diversity table
meta_div.to_csv("diversity_indices.tsv", sep="\t", index=False)
print("Diversity indices saved to diversity_indices.tsv")

# --- Visualization ---

sns.set(style="whitegrid", font_scale=1.2)

# Boxplot: Shannon vs study_condition
plt.figure(figsize=(6,5))
sns.boxplot(x="study_condition", y="Shannon", data=meta_div, palette="Set2")
sns.stripplot(x="study_condition", y="Shannon", data=meta_div, color="black", alpha=0.5)
plt.title("Shannon Index by Study Condition")
plt.tight_layout()
plt.savefig("shannon_by_condition.png", dpi=300)

# Boxplot: Pielou vs study_condition
plt.figure(figsize=(6,5))
sns.boxplot(x="study_condition", y="Pielou", data=meta_div, palette="Set2")
sns.stripplot(x="study_condition", y="Pielou", data=meta_div, color="black", alpha=0.5)
plt.title("Pielou's Evenness by Study Condition")
plt.tight_layout()
plt.savefig("pielou_by_condition.png", dpi=300)

# Boxplot: Shannon vs age_category
plt.figure(figsize=(7,5))
sns.boxplot(x="age_category", y="Shannon", data=meta_div, palette="coolwarm")
sns.stripplot(x="age_category", y="Shannon", data=meta_div, color="black", alpha=0.5)
plt.title("Shannon Index by Age Category")
plt.tight_layout()
plt.savefig("shannon_by_age.png", dpi=300)

# Boxplot: Pielou vs age_category
plt.figure(figsize=(7,5))
sns.boxplot(x="age_category", y="Pielou", data=meta_div, palette="coolwarm")
sns.stripplot(x="age_category", y="Pielou", data=meta_div, color="black", alpha=0.5)
plt.title("Pielou's Evenness by Age Category")
plt.tight_layout()
plt.savefig("pielou_by_age.png", dpi=300)

print("All Q-4 boxplots saved.")
EOF
