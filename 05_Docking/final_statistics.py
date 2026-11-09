import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# 1. Load Data
# Ensure this CSV contains ALL 33 ligands and their best affinities
# If you don't have a single CSV with all best scores, you can create one manually 
# or use the 'docking_results_summary.csv' if it has all entries.
csv_path = os.path.join(os.path.dirname(__file__), 'docking_results_summary.csv')
df = pd.read_csv(csv_path)

# Clean Ligand Names (Remove conformer info like _0_1)
df['Ligand_Name'] = df['Filename'].apply(lambda x: x.split('_0_')[0].replace('.pdbqt', ''))

# Get the BEST (lowest) affinity for each unique ligand
best_scores = df.groupby('Ligand_Name')['Affinity_kcal_mol'].min().reset_index()

# 2. Define Groups
top_hits_list = ['Apigenin', 'Kaempferol', 'Luteolin', 'Quercetin', 'Scutellarein']
reference_list = ['RCX']

# Extract Data for Stats
hit_data = best_scores[best_scores['Ligand_Name'].isin(top_hits_list)]['Affinity_kcal_mol']
ref_data = best_scores[best_scores['Ligand_Name'].isin(reference_list)]['Affinity_kcal_mol']
rest_data = best_scores[~best_scores['Ligand_Name'].isin(top_hits_list + reference_list)]['Affinity_kcal_mol']

# 3. Statistical Tests
print("--- Statistical Analysis ---")

# A. T-Test: Are Top Hits significantly better (lower energy) than the Rest?
t_stat, p_value_hits_vs_rest = stats.ttest_ind(hit_data, rest_data)
print(f"1. T-Test (Top Hits vs. Rest of Library):")
print(f"   Mean Hit Affinity: {hit_data.mean():.2f} kcal/mol")
print(f"   Mean Rest Affinity: {rest_data.mean():.2f} kcal/mol")
print(f"   P-Value: {p_value_hits_vs_rest:.4e}")
if p_value_hits_vs_rest < 0.05:
    print("   Result: SIGNIFICANT (Hits are statistically better than the rest)")
else:
    print("   Result: Not Significant")

# B. Comparison to Reference (RCX)
ref_val = ref_data.values[0]
print(f"\n2. Comparison to Reference (RCX: {ref_val:.2f} kcal/mol):")
for name in top_hits_list:
    score = best_scores[best_scores['Ligand_Name'] == name]['Affinity_kcal_mol'].values[0]
    diff = score - ref_val
    print(f"   {name}: {score:.2f} kcal/mol (Diff: {diff:+.2f})")

# 4. Visualization: Histogram Plot
mean_affinity = best_scores['Affinity_kcal_mol'].mean()
better_count = (best_scores['Affinity_kcal_mol'] < ref_val).sum()

plt.figure(figsize=(12, 8))
ax = sns.histplot(
    best_scores['Affinity_kcal_mol'],
    bins=15,
    kde=False,
    color='#6baed6',
    edgecolor='black',
    linewidth=0.8,
)

# Reference and mean lines
ref_line = ax.axvline(ref_val, color='green', linestyle='--', linewidth=2, label=f'Reference Ligand (RCX): {ref_val:.2f}')
mean_line = ax.axvline(mean_affinity, color='red', linestyle=':', linewidth=2, label=f'Average Binding Affinity = {mean_affinity:.2f}')

# Shade compounds better than reference
x_min, x_max = ax.get_xlim()
ax.axvspan(x_min, ref_val, color='green', alpha=0.12)

# Legend patch for shaded region
import matplotlib.patches as mpatches
shade_patch = mpatches.Patch(color='green', alpha=0.12, label='Compounds Better than Reference')
ax.legend(handles=[shade_patch, ref_line, mean_line], loc='upper left')

# Labels and styling
ax.set_title('Binding Affinity Distribution of Dicoma anomala Phytochemicals', fontsize=16)
ax.set_xlabel('Binding Affinity (kcal/mol)', fontsize=14)
ax.set_ylabel('Number of Ligands', fontsize=14)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

# Tight layout and save
plt.tight_layout()
plt.savefig('binding_affinity_distribution.png', dpi=300)
plt.show()

print(f"\nPlot saved as 'binding_affinity_distribution.png' ({better_count} compounds better than RCX)")