#!/bin/bash
# Q-3 SPINGO pipeline for Batch3

# Paths
SPINGO_PATH=~/Batch3_16S/SPINGO-master/spingo
RDP_DB=~/Batch3_16S/SPINGO-master/database/RDP_database/database/RDP_11.2.species.fa
THREADS=10

# Directories
mkdir -p spingo_fasta spingo_out

echo "Step 1: Running SPINGO on all paired samples..."

# Loop through all paired trimmed samples
for f1 in trimmed/*_1_paired.fastq.gz; do
    base=$(basename $f1 _1_paired.fastq.gz)
    f2=trimmed/${base}_2_paired.fastq.gz

    echo "Processing sample: $base"

    # concatenate paired reads
    zcat $f1 $f2 > spingo_fasta/${base}.fastq

    # convert fastq -> fasta
    awk 'NR%4==1 {gsub("@",">",$0); print} NR%4==2 {print}' spingo_fasta/${base}.fastq > spingo_fasta/${base}.fasta

    # run SPINGO
    $SPINGO_PATH -d $RDP_DB -p $THREADS -i spingo_fasta/${base}.fasta > spingo_out/${base}.spingo.out.txt

    # clean intermediates
    rm spingo_fasta/${base}.fastq spingo_fasta/${base}.fasta
done

echo "Step 2: Parsing SPINGO outputs to create counts table..."

# Python inline for parsing and plotting
python3 << 'EOF'
import glob, pandas as pd
import matplotlib.pyplot as plt

files = sorted(glob.glob('spingo_out/*.spingo.out.txt'))
all_species = {}

for fn in files:
    sample = fn.split('/')[-1].replace('.spingo.out.txt','')
    counts = {}
    with open(fn) as fh:
        for line in fh:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            species, species_conf = cols[-2], float(cols[-1])
            # filter by confidence >= 0.7 and remove ambiguous
            if species_conf >= 0.7 and species.lower() != 'ambiguous':
                counts[species] = counts.get(species,0)+1
    all_species[sample] = counts

# Create counts dataframe
species_list = sorted({s for d in all_species.values() for s in d.keys()})
df = pd.DataFrame(0, index=species_list, columns=all_species.keys())
for samp, d in all_species.items():
    for sp, c in d.items():
        df.at[sp, samp] = c

df.to_csv('species_counts_table.tsv', sep='\t')
print("species_counts_table.tsv created.")

# Total-sum scaling normalization
rel_ab = df.div(df.sum(axis=0), axis=1)
rel_ab.to_csv('species_rel_abundance.tsv', sep='\t')
print("species_rel_abundance.tsv created.")

# Top 20 species by mean relative abundance
top20 = rel_ab.mean(axis=1).sort_values(ascending=False).head(20)
top20.to_csv('top20_species_mean_abundance.tsv', sep='\t', header=['mean_abundance'])

# Plot top 20 species
top20.sort_values().plot(kind='barh', figsize=(6,8), color='skyblue')
plt.xlabel('Mean Relative Abundance')
plt.ylabel('Species')
plt.title('Top 20 Species (Batch3)')
plt.tight_layout()
plt.savefig('top20_species_barplot.png', dpi=300)
# plt.show()
EOF

echo "Q-3 pipeline finished. Outputs:"
echo "- species_counts_table.tsv"
echo "- species_rel_abundance.tsv"
echo "- top20_species_mean_abundance.tsv"
echo "- top20_species_barplot.png"
