# 16S rRNA Microbiome Analysis Pipeline - Batch3

## Overview
This pipeline processes 16S rRNA sequencing data through quality control, trimming, taxonomic classification, and diversity analysis for 20 samples (40 FASTQ files). The analysis includes quality assessment, adapter trimming, species-level classification using SPINGO, and diversity analysis with Shannon and Pielou indices.

## Dependencies
- **FastQC** v0.11.9
- **Trimmomatic** v0.39
- **SPINGO** v1.3
- **Python** 3.8+ with packages: pandas, matplotlib, seaborn, numpy
- **RDP database** v11.2

## Project Structure
```
Batch3_16S/
├── fastqc_raw/                 # Raw read quality reports
├── trimmed/                    # Trimmed reads
├── fastqc_trimmed/            # Post-trimming quality reports
├── spingo_fasta/              # Temporary FASTA files
├── spingo_out/                # SPINGO classification results
├── Batch3_samples.sh          # Sample download script
├── Batch3.csv                 # Metadata file
├── run_spingo_Q3.sh           # SPINGO analysis script
└── run_Q4_diversity.sh        # Diversity analysis script
```

## Pipeline Steps

### Question 1: Quality Control of Raw Reads

```bash
# Create working directory
mkdir -p ~/Batch3_16S && cd ~/Batch3_16S

# Copy necessary files
cp ~/Downloads/Batch3_samples.sh .
cp ~/Downloads/Batch3.csv .

# Download samples and run FastQC
bash Batch3_samples.sh

# Activate FastQC environment and run quality control
conda activate fastqc_env
mkdir fastqc_raw
fastqc -t 4 -o fastqc_raw *.fastq.gz
multiqc -o fastqc_raw fastqc_raw
```

**Outputs:**
- Raw read quality reports: `fastqc_raw/`
- Combined quality summary: `fastqc_raw/multiqc_report.html`

### Question 2: Trimming and Post-Trim QC

```bash
Step-1:
mkdir -p trimmed
conda activate trimmomatic_env

Step-2:
# Setting path to adapter file
ADAPTERS=$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE.fa
for f1 in *_1.fastq.gz; do
f2=${f1/_1.fastq.gz/_2.fastq.gz}
base=${f1%_1.fastq.gz}
trimmomatic PE -threads 4 -phred33 \
$f1 $f2 \
trimmed/${base}_1_paired.fastq.gz trimmed/${base}_1_unpaired.fastq.gz \
trimmed/${base}_2_paired.fastq.gz trimmed/${base}_2_unpaired.fastq.gz \
ILLUMINACLIP:${ADAPTERS}:2:30:10 \
SLIDINGWINDOW:5:27 AVGQUAL:27 MINLEN:100
done

Step-3:
# Post-trimming quality check
mkdir fastqc_trimmed
fastqc -t 4 -o fastqc_trimmed trimmed/*_paired.fastq.gz
multiqc -o fastqc_trimmed fastqc_trimmed
# Post-trimming quality check
mkdir fastqc_trimmed
fastqc -t 4 -o fastqc_trimmed trimmed/*_paired.fastq.gz
multiqc -o fastqc_trimmed fastqc_trimmed
```

**Outputs:**
- High-quality paired reads: `trimmed/*_paired.fastq.gz`
- Post-trim quality reports: `fastqc_trimmed/`

### Question 3: Taxonomic Classification with SPINGO

```bash
#!/bin/bash
# Q-3 SPINGO pipeline for Batch3
# This script automates species-level classification and summarization for all samples in Batch3.

# --------------------------------------------------------
# Define important paths and parameters
# --------------------------------------------------------
SPINGO_PATH=~/Batch3_16S/SPINGO-master/spingo             # Path to the SPINGO executable
RDP_DB=~/Batch3_16S/SPINGO-master/database/RDP_database/database/RDP_11.2.species.fa   # Path to the RDP reference database used for classification
THREADS=10                                                # Number of CPU threads to use (I have 12 total; using 10 leaving some free)

# --------------------------------------------------------
# Create required directories for intermediate and output files
# --------------------------------------------------------
mkdir -p spingo_fasta spingo_out                          # Create folders to store temporary FASTA files and SPINGO output results

echo "Step 1: Running SPINGO on all paired samples..."

# --------------------------------------------------------
# Loop through all paired-end trimmed FASTQ files
# --------------------------------------------------------
for f1 in trimmed/*_1_paired.fastq.gz; do                 # Iterate through each forward (R1) trimmed FASTQ file
    base=$(basename $f1 _1_paired.fastq.gz)               # Extract the base sample name (e.g., ERR4042824)
    f2=trimmed/${base}_2_paired.fastq.gz                  # Define the corresponding reverse (R2) file path

    echo "Processing sample: $base"

    # ----------------------------------------------------
    # Step 1a: Concatenate forward and reverse reads
    # ----------------------------------------------------
    zcat $f1 $f2 > spingo_fasta/${base}.fastq             # Combine both R1 and R2 into a single FASTQ file (SPINGO expects a single input file)

    # ----------------------------------------------------
    # Step 1b: Convert FASTQ → FASTA (SPINGO requires FASTA input)
    # ----------------------------------------------------
    awk 'NR%4==1 {gsub("@",">",$0); print} NR%4==2 {print}' \
        spingo_fasta/${base}.fastq > spingo_fasta/${base}.fasta

    # ----------------------------------------------------
    # Step 1c: Run SPINGO classification
    # ----------------------------------------------------
    $SPINGO_PATH -d $RDP_DB -p $THREADS -i spingo_fasta/${base}.fasta \
        > spingo_out/${base}.spingo.out.txt               # Run SPINGO using the RDP database and save results per sample

    # ----------------------------------------------------
    # Step 1d: Clean up temporary intermediate files
    # ----------------------------------------------------
    rm spingo_fasta/${base}.fastq spingo_fasta/${base}.fasta  # Delete temporary FASTQ/FASTA files to save space
done

echo "Step 2: Parsing SPINGO outputs to create counts table..."

# --------------------------------------------------------
# Step 3: Inline Python script for parsing & summarization
# --------------------------------------------------------
python3 << 'EOF'
import glob, pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------------
# Collect all SPINGO result files
# --------------------------------------------------------
files = sorted(glob.glob('spingo_out/*.spingo.out.txt'))  # Locate all per-sample SPINGO output files
all_species = {}                                          # Dictionary to store species counts per sample

# --------------------------------------------------------
# Parse each SPINGO output
# --------------------------------------------------------
for fn in files:
    sample = fn.split('/')[-1].replace('.spingo.out.txt','')   # Extract sample name
    counts = {}
    with open(fn) as fh:
        for line in fh:
            if line.startswith('#'): continue                  # Skip comment lines
            cols = line.strip().split('\t')                    # Split columns by tab
            species, species_conf = cols[-2], float(cols[-1])  # Extract species name and confidence score
            # Filter by confidence >= 0.7 and exclude 'ambiguous' classifications
            if species_conf >= 0.7 and species.lower() != 'ambiguous':
                counts[species] = counts.get(species,0)+1
    all_species[sample] = counts                               # Store counts for this sample

# --------------------------------------------------------
# Create a full species count table
# --------------------------------------------------------
species_list = sorted({s for d in all_species.values() for s in d.keys()})  # Unique list of all detected species
df = pd.DataFrame(0, index=species_list, columns=all_species.keys())        # Initialize dataframe with zeros
for samp, d in all_species.items():                                         # Fill table with counts
    for sp, c in d.items():
        df.at[sp, samp] = c

df.to_csv('species_counts_table.tsv', sep='\t')                            # Save raw counts
print("species_counts_table.tsv created.")

# --------------------------------------------------------
# Total-sum scaling (TSS) normalization
# --------------------------------------------------------
rel_ab = df.div(df.sum(axis=0), axis=1)                                   # Normalize each column by its total count
rel_ab.to_csv('species_rel_abundance.tsv', sep='\t')                      # Save relative abundance table
print("species_rel_abundance.tsv created.")

# --------------------------------------------------------
# Identify top 20 most abundant species
# --------------------------------------------------------
top20 = rel_ab.mean(axis=1).sort_values(ascending=False).head(20)         # Compute mean relative abundance across samples
top20.to_csv('top20_species_mean_abundance.tsv', sep='\t', header=['mean_abundance'])  # Save top 20 list

# --------------------------------------------------------
# Plot top 20 species bar chart
# --------------------------------------------------------
top20.sort_values().plot(kind='barh', figsize=(6,8), color='skyblue')     # Create horizontal bar plot
plt.xlabel('Mean Relative Abundance')
plt.ylabel('Species')
plt.title('Top 20 Species (Batch3)')
plt.tight_layout()
plt.savefig('top20_species_barplot.png', dpi=300)                         # Save plot as PNG
# plt.show()  # Disabled to prevent blocking execution in non-GUI mode
EOF

# --------------------------------------------------------
# Step 4: Print final output summary
# --------------------------------------------------------
echo "Q-3 pipeline finished. Outputs:"
echo "- species_counts_table.tsv"
echo "- species_rel_abundance.tsv"
echo "- top20_species_mean_abundance.tsv"
echo "- top20_species_barplot.png"


## in Trimmomatric Environment
chmod +x ./run_spingo_Q3.sh
./run_spingo_Q3.sh
```

**Key Script Components:**
- Concatenates forward and reverse reads
- Converts FASTQ to FASTA format
- Runs SPINGO classification with RDP database
- Generates species abundance tables and visualizations

**Outputs:**
- Raw counts: `species_counts_table.tsv`
- Relative abundances: `species_rel_abundance.tsv`
- Top 20 species: `top20_species_mean_abundance.tsv`
- Visualization: `top20_species_barplot.png`

### Question 4: Diversity Analysis

```bash
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


# Running the Code
chmod +x ./run_Q4_diversity.sh && ./run_Q4_diversity.sh
```

**Analysis Includes:**
- Shannon diversity index calculation
- Pielou's evenness measurement
- Statistical comparisons by study condition and age category
- Visualization of diversity patterns

**Outputs:**
- Diversity indices: `diversity_indices.tsv`
- Boxplots:
  - `shannon_by_condition.png`
  - `pielou_by_condition.png`
  - `shannon_by_age.png`
  - `pielou_by_age.png`

## Results Summary

The pipeline successfully processes 20 paired-end samples through:
1. **Quality Control**: Identifies and removes low-quality sequences
2. **Trimming**: Removes adapters and improves read quality
3. **Taxonomic Classification**: Assigns species-level taxonomy using RDP database
4. **Diversity Analysis**: Computes alpha diversity metrics and visualizes patterns

All intermediate files and final results are organized in dedicated directories for reproducibility and easy access.

## Group Information
- **Group**: 3
- **Date**: 7th October 2024
