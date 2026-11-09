import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Load Data
csv_path = os.path.join(os.path.dirname(__file__), 'docking_results_summary.csv')
df = pd.read_csv(csv_path)

# Clean Ligand Names
df['Ligand_Name'] = df['Filename'].apply(lambda x: x.split('_0_')[0].replace('.pdbqt', ''))

# Get best scores
best_scores = df.groupby('Ligand_Name')['Affinity_kcal_mol'].min().reset_index()

# Define Groups
top_hits_list = ['Apigenin', 'Kaempferol', 'Luteolin', 'Quercetin', 'Scutellarein']
reference_list = ['RCX']

# Extract Data
hit_data = best_scores[best_scores['Ligand_Name'].isin(top_hits_list)]['Affinity_kcal_mol']
ref_data = best_scores[best_scores['Ligand_Name'].isin(reference_list)]['Affinity_kcal_mol']
rest_data = best_scores[~best_scores['Ligand_Name'].isin(top_hits_list + reference_list)]['Affinity_kcal_mol']

print("--- Advanced Statistical Analysis ---")

# Descriptive Statistics
print("\n1. Descriptive Statistics:")
for group_name, data in [("Top Hits", hit_data), ("Reference (RCX)", ref_data), ("Rest of Library", rest_data)]:
    print(f"\n{group_name}:")
    print(f"   Count: {len(data)}")
    print(f"   Mean: {data.mean():.2f} kcal/mol")
    print(f"   Std Dev: {data.std():.2f}")
    print(f"   Median: {data.median():.2f}")
    print(f"   Min: {data.min():.2f}")
    print(f"   Max: {data.max():.2f}")

# ANOVA: One-way ANOVA across all groups
all_groups = [hit_data, ref_data, rest_data]
f_stat, p_value_anova = stats.f_oneway(*all_groups)
print(f"\n2. One-Way ANOVA (Top Hits vs. Reference vs. Rest):")
print(f"   F-Statistic: {f_stat:.2f}")
print(f"   P-Value: {p_value_anova:.4e}")
if p_value_anova < 0.05:
    print("   Result: SIGNIFICANT differences between groups")
else:
    print("   Result: No significant differences between groups")

# Post-hoc: Tukey's HSD if ANOVA significant
if p_value_anova < 0.05:
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    all_data = pd.concat([hit_data, ref_data, rest_data])
    labels = ['Top Hits'] * len(hit_data) + ['Reference'] * len(ref_data) + ['Rest'] * len(rest_data)
    tukey = pairwise_tukeyhsd(endog=all_data, groups=labels, alpha=0.05)
    print(f"\n3. Post-hoc Tukey's HSD:")
    print(tukey)

# Correlation Analysis
# Since we only have affinity data, we'll compute correlation matrix for all affinities
# But with one variable, correlation is 1. Perhaps correlate with rank or something, but let's note.
print(f"\n4. Correlation Analysis:")
print("   Note: With only binding affinity data available, full correlation analysis requires additional variables (e.g., molecular weight, logP).")
print("   Correlation matrix for affinities (trivial, all 1.0):")
corr_matrix = best_scores[['Affinity_kcal_mol']].corr()
print(corr_matrix)

# Perhaps correlate affinity with the number of conformers per ligand (if available)
# But from data, we can count how many entries per ligand
conformer_count = df.groupby('Ligand_Name').size().reset_index(name='Conformer_Count')
best_scores_with_count = best_scores.merge(conformer_count, on='Ligand_Name')
corr_affinity_conformers = best_scores_with_count[['Affinity_kcal_mol', 'Conformer_Count']].corr()
print(f"\n   Correlation between Affinity and Number of Conformers:")
print(corr_affinity_conformers)

print("\nAnalysis complete.")