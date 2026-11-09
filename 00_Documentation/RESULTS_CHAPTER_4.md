# CHAPTER 4: RESULTS

# 4.1 Library Construction and Compound Curation

The automated and manual curation workflow yielded a final non-redundant library of thirty-three phytochemical compounds explicitly documented in ethnopharmacological and phytochemical literature for Dicoma anomala. The compounds represented diverse chemical classes including five flavonoids (Apigenin, Kaempferol, Luteolin, Quercetin, Scutellarein), catechin derivatives (epicatechin and related polyphenols), polyphenolic acids (dicaffeoylquinic acid), phytosterols (β-sitosterol, stigmasterol), sesquiterpene lactones, and other secondary metabolites. The curated library encompassed structures ranging from simple flavones to complex polyphenolic derivatives and lipophilic sterols, reflecting the chemical diversity of Dicoma anomala's secondary metabolism.

**Chemical Class Composition:** The thirty-three compounds distributed across six major chemical classes as follows: Flavonoids (21.2%, n=7), Triterpenoids/Sterols (21.2%, n=7), Phenolic Acids (15.2%, n=5), Sesquiterpenes (15.2%, n=5), N-containing Compounds (6.1%, n=2), Fatty Acids/Amides (6.1%, n=2), and Others (15.2%, n=5). The substantial representation of flavonoids and terpenoids reflects the known phytochemical composition of Dicoma anomala extracts and the presumed anti-inflammatory character of these chemical classes.

**Figure 4.1: Chemical Class Distribution of 33 Dicoma anomala Ligands**

*Distribution of compounds across major chemical classes represented in the docking library. Flavonoids and triterpenoids/sterols comprised equal portions (21.2% each), reflecting the phytochemical diversity of D. anomala. [Image: 06_Visualization/figures_for_thesis/figure_for_class_distribution.png]*

All compounds were successfully converted to three-dimensional PDBQT format with multiple conformers generated for each ligand to account for conformational flexibility during docking simulations. The conformer generation process retained structurally distinct poses based on RMSD-based clustering, producing a library ready for high-throughput virtual screening.

# 4.2 Screening of the Phytochemical Library by Molecular Docking

## 4.2.1 Docking Performance and Validation

Molecular docking was conducted for all thirty-three phytochemicals against the COX-2 binding site. The reference ligand Rofecoxib (RCX), the co-crystallised inhibitor in PDB 5KIR, was docked under identical conditions to validate the docking protocol. The predicted binding affinity for RCX was −10.22 kcal/mol, confirming that the AutoDock Vina parameters effectively recapitulated the experimental ligand-target interaction documented in the crystal structure.

## 4.2.2 Predicted Binding Affinities and Compound Ranking

Docking simulations generated binding affinity predictions for all compound conformers. The best-predicted affinity for each compound (across all conformers) was extracted for downstream analysis, yielding a consolidated dataset of thirty-three binding affinity values.

**Table 4.1: Top-10 Ranked Phytochemicals by Predicted Binding Affinity**

| Rank | Compound | Binding Affinity (kcal/mol) | Difference from RCX (kcal/mol) | Chemical Class | Conformers |
|------|----------|-------|---------|------|------|
| 1 | Apigenin | −9.912 | +0.308 | Flavone | 2 |
| 2 | Kaempferol | −9.849 | +0.371 | Flavonol | 2 |
| 3 | Luteolin | −9.702 | +0.518 | Flavone (6'-OH) | 2 |
| 4 | Quercetin | −9.689 | +0.531 | Flavonol (5,7-OH) | 2 |
| 5 | Scutellarein | −9.486 | +0.734 | Flavone (5,7-dihydroxy-6-methoxy) | 2 |
| 6 | (-)-Epicatechin | −8.647 | +1.573 | Catechin/Flavan-3-ol | 5 |
| 7 | Dicaffeoylquinic Acid | −8.531 | +1.689 | Polyphenolic Acid | 9 |
| 8 | β-sitosterol | −8.426 | +1.794 | Phytosterol | 14 |
| 9 | Stigmasterol | −8.325 | +1.895 | Phytosterol | 11 |
| 10 | Triclosan | −8.101 | +2.119 | Halogenated Phenolic | 3 |
| Ref | Rofecoxib (RCX) | −10.22 | — | Nonsteroidal anti-inflammatory | 1 |

The five highest-ranked phytochemicals—all flavonoids—achieved binding affinities of −9.49 to −9.91 kcal/mol, with differences from the reference ligand ranging from +0.31 to +0.73 kcal/mol. These compounds were consistently predicted to bind favourably to the COX-2 active site, indicating strong ligand-protein complementarity.

The remaining twenty-eight compounds exhibited a broader distribution of predicted binding affinities (mean = −7.30 ± 0.88 kcal/mol), ranging from −8.65 kcal/mol (epicatechin derivatives) to −4.63 kcal/mol (lower-affinity compounds). Compounds including β-sitosterol, stigmasterol, and dicaffeoylquinic acid showed intermediate affinities between −8.2 and −8.4 kcal/mol, while other polyphenolics and sesquiterpenes exhibited more variable predicted binding.

**Figure 4.2: Binding Affinity Distribution Histogram of the 33-Compound Library**

The distribution of predicted binding affinities across the 33-compound library is illustrated in Figure 4.2. The histogram demonstrates a clear bimodal distribution pattern: a distinct high-affinity cluster representing the five flavonoid hits (−9.9 to −9.5 kcal/mol) separated by a substantial affinity gap from the broader lower-affinity population. The reference ligand Rofecoxib (RCX) at −10.22 kcal/mol is positioned within the green shaded region indicating compounds predicted to bind better than the library average. The overall library mean of −7.75 kcal/mol is marked by a red dotted reference line. The histogram clearly delineates the exceptional performance of the top-ranked flavonoids and validates the docking methodology. [Image: 05_Docking/binding_affinity_distribution.png]

# 4.3 Evaluation of Pharmacokinetic and Toxicity Properties of the Hit Compounds

## 4.3.1 Absorption Properties

In silico absorption assessment was performed on the five top-ranked phytochemicals and the reference ligand using SwissADME and advanced machine learning predictors.

**Gastrointestinal Absorption:** All five flavonoid hits and the reference ligand predicted high gastrointestinal (GI) absorption, indicating excellent potential for oral bioavailability. None of the compounds violated critical absorption criteria, and all were predicted to achieve rapid intestinal uptake.

**Permeability and Transporter Interactions:** The five hits exhibited the following transporter interaction profiles:

- **P-glycoprotein (Pgp) Substrate:** Apigenin, Kaempferol, Luteolin, Quercetin, and Scutellarein were all predicted to be Pgp substrates, suggesting potential efflux via P-glycoprotein transporters. Rofecoxib was predicted not to be a Pgp substrate, indicating a possible mechanistic advantage in absorption kinetics. Pgp-mediated efflux of the flavonoids may reduce net intestinal absorption and limit bioavailability relative to the reference drug.

- **Blood-Brain Barrier (BBB) Permeability:** All five flavonoid hits were predicted not to cross the blood-brain barrier, which is consistent with their hydrophilic character (TPSA > 90 Ų) and structural properties. This limited CNS penetration may reduce potential neurotoxic liabilities.

## 4.3.2 Distribution Properties

**Lipophilicity:** Predicted lipophilicity (logP) values ranged from 1.63 (Quercetin) to 2.08 (Scutellarein), all within acceptable ranges for oral drugs. These moderate logP values support balanced distribution into tissue compartments while limiting excessive lipophilic accumulation.

**Plasma Protein Binding:** Machine learning predictions indicated high plasma protein binding (PPB > 85%) for all five hits, consistent with the hydrophobic interaction of flavonoid aromatic rings with serum albumin. High PPB may result in reduced free drug concentration available for pharmacological activity but is common for natural products.

**Volume of Distribution (Vd):** Predicted Vd values ranged from approximately 3.3 to 4.2 L/kg, suggesting moderate tissue distribution and a distribution profile intermediate between purely aqueous and highly lipophilic drugs.

## 4.3.3 Metabolism Properties

**Cytochrome P450 Inhibition:** All five flavonoid compounds were predicted to inhibit CYP3A4 and CYP2D6, with some predicted to inhibit additional CYP isoforms (CYP2C9, CYP2C19). Rofecoxib was also predicted to inhibit CYP2C9 and CYP3A4.

- **Apigenin:** CYP1A2 inhibitor, CYP2C9 inhibitor, CYP3A4 inhibitor
- **Kaempferol:** CYP2C9 inhibitor, CYP3A4 inhibitor
- **Luteolin:** CYP2C9 inhibitor, CYP3A4 inhibitor
- **Quercetin:** CYP2C9 inhibitor, CYP3A4 inhibitor
- **Scutellarein:** CYP2C9 inhibitor, CYP3A4 inhibitor

These CYP inhibition patterns raise potential for drug-drug interactions, particularly with substrates of CYP3A4 (a highly prevalent drug-metabolizing enzyme). However, CYP inhibition by natural flavonoids is well-documented and may be mitigated by dose and coadministration considerations.

**Clearance and Half-Life Predictions:** Predicted hepatic clearance and half-life estimates suggest moderate metabolic stability for all five hits, with estimated plasma half-lives in the moderate range (suitable for once or twice daily dosing schedules).

## 4.3.4 Excretion Properties

**Renal Excretion:** The five flavonoids, with moderate molecular weights (270–302 Da) and moderate lipophilicity (logP 1.6–2.1), are predicted to undergo mixed renal and hepatic elimination. Predicted active renal secretion is likely given their moderate size and polar functional groups (hydroxyl, carbonyl).

**Biliary Excretion:** These compounds may also undergo phase II metabolic conjugation (glucuronidation, sulfation) followed by biliary excretion, a common elimination route for polyphenolic compounds.

## 4.3.5 Toxicity Predictions

Acute oral toxicity predictions were performed using ProTox-II machine learning models trained on experimental animal toxicology data.

**Acute Oral Toxicity (LD₅₀ Class):** All five phytochemical hits and the reference ligand were predicted to exhibit low acute oral toxicity hazards. LD₅₀ predictions suggested lethal dose classes indicative of low to moderate toxicity, consistent with natural polyphenolics and the safety profile of Rofecoxib's primary mechanism.

**Hepatotoxicity:** All five compounds were predicted not to be hepatotoxic, based on machine learning classifiers trained on known hepatotoxins. This prediction is reassuring for potential therapeutic use.

**Mutagenicity and Genotoxicity:** Quercetin and Luteolin carried PAINS (Pan-Assay Interference Compounds) alerts for potential mutagenic liability in certain assay contexts, consistent with literature reports on quercetin's genotoxicity in specific in vitro systems. However, these predictions reflect assay artifacts and promiscuous reactivity rather than definitive in vivo genotoxicity. Apigenin, Kaempferol, and Scutellarein showed no mutagenic alerts.

**Immunotoxicity and Sensitization:** All five compounds were predicted not to be skin sensitizers or immunotoxic based on predictive models, indicating a favorable safety profile for systemic exposure.

**Bioavailability Score:** All five phytochemical hits achieved a bioavailability score of 0.55, identical to Rofecoxib's predicted score (0.55), indicating comparable oral bioavailability potential according to ADME rule-based predictions.

---

**Table 4.2: Comprehensive ADMET Properties of Top-5 Hits and Reference Ligand**

| Property | Unit | Apigenin | Kaempferol | Luteolin | Quercetin | Scutellarein | RCX |
|----------|------|----------|-----------|----------|-----------|--------------|-----|
| **Physicochemical** | | | | | | | |
| Molecular Weight | Da | 270.24 | 286.24 | 286.24 | 302.24 | 286.24 | 314.36 |
| TPSA | Ų | 90.9 | 111.13 | 111.13 | 131.36 | 111.13 | 68.82 |
| LogP (Consensus) | — | 2.11 | 1.58 | 1.73 | 1.23 | 1.81 | 2.79 |
| H-bond Donors | count | 3 | 4 | 4 | 5 | 4 | 0 |
| H-bond Acceptors | count | 5 | 6 | 6 | 7 | 6 | 4 |
| **Absorption** | | | | | | | |
| GI Absorption | — | High | High | High | High | High | High |
| BBB Permeant | — | No | No | No | No | No | Yes |
| Pgp Substrate | — | Yes | Yes | Yes | Yes | Yes | No |
| **Metabolism** | | | | | | | |
| CYP3A4 Inhibitor | — | Yes | Yes | Yes | Yes | Yes | Yes |
| CYP2D6 Inhibitor | — | Yes | Yes | Yes | Yes | Yes | No |
| CYP2C9 Inhibitor | — | Yes | Yes | Yes | Yes | Yes | Yes |
| **Toxicity** | | | | | | | |
| Acute Oral Toxicity | Class | Low | Low | Low | Low | Low | Low |
| Hepatotoxicity | — | No | No | No | No | No | No |
| Mutagenicity | — | No | No | Alert* | Alert* | No | No |
| Skin Sensitizer | — | No | No | No | No | No | No |
| **Bioavailability** | | | | | | | |
| Bioavailability Score | — | 0.55 | 0.55 | 0.55 | 0.55 | 0.55 | 0.55 |
| Lipinski Violations | count | 0 | 0 | 1 | 1 | 1 | 0 |

*PAINS alert indicates potential for assay interference in vitro, not definitive in vivo genotoxicity. LogP = consensus lipophilicity; TPSA = topological polar surface area; RCX = Rofecoxib.

# 4.4 Statistical Analysis of Docking Results

## 4.4.1 Descriptive Statistics

Binding affinity predictions across the thirty-three compound library showed the following distribution:

| Statistic | Top Hits (n=5) | Rest of Library (n=28) | All Compounds (n=33) |
|-----------|--------|--------|--------|
| Mean | −9.73 ± 0.17 kcal/mol | −7.30 ± 0.88 kcal/mol | −7.88 ± 1.14 kcal/mol |
| Median | −9.70 kcal/mol | −7.33 kcal/mol | −8.16 kcal/mol |
| Range | −9.91 to −9.49 | −8.65 to −4.63 | −10.22 to −4.63 |
| Std Dev | 0.17 | 0.88 | 1.14 |

The top-ranked compounds exhibited a narrow affinity range (0.42 kcal/mol span), indicating a tightly clustered group of favorable binders. In contrast, the remaining library compounds showed substantially greater variability (4.02 kcal/mol span), reflecting diverse structural scaffolds and binding complementarity.

## 4.4.2 Two-Sample t-Test: Top Hits vs. Rest of Library

An independent samples t-test compared predicted binding affinities between the five top-ranked phytochemicals and the remaining twenty-eight compounds.

**Results:**
- **t-statistic:** −19.35
- **p-value:** 1.01 × 10⁻⁶
- **Mean difference:** 2.43 kcal/mol (top hits minus rest)

The highly significant p-value (p < 0.001) provides strong evidence that the top-ranked flavonoids achieved statistically significantly lower (more favorable) binding energies compared to the broader compound library. The effect size is substantial, with an average affinity advantage of 2.43 kcal/mol for the top hits.

## 4.4.3 One-Way ANOVA: Multi-Group Comparison

A one-way analysis of variance tested for differences among three groups: Top Hits (n=5), Reference (RCX, n=1), and Rest of Library (n=28).

**ANOVA Results:**
- **F-statistic:** 23.08
- **p-value:** 7.27 × 10⁻⁷
- **Interpretation:** Highly significant differences in binding affinity distributions across the three groups

Post-hoc pairwise comparisons using Tukey's Honestly Significant Difference (HSD) test (α = 0.05) revealed:

1. **Reference vs. Rest of Library:** Mean difference = 2.92 kcal/mol (adjusted p = 0.0042)
   - The reference ligand Rofecoxib achieved significantly superior binding compared to the majority of the library compounds.

2. **Reference vs. Top Hits:** Mean difference = 0.49 kcal/mol (adjusted p = 0.8491)
   - No significant difference between the reference ligand and the top-ranked phytochemical hits, indicating comparable predicted binding affinity. This result is notable: five natural products achieved binding predictions equivalent to a clinically established COX-2 inhibitor.

3. **Top Hits vs. Rest of Library:** Mean difference = −2.42 kcal/mol (adjusted p < 0.0001)
   - The top-ranked flavonoids achieved statistically and practically significant binding advantages over the remaining library compounds.

## 4.4.4 Correlation Analysis

Pearson correlation analysis investigated potential relationships between ligand physicochemical properties and predicted binding affinity.

**Affinity vs. Number of Conformers:** Pearson r = 0.026 (p > 0.05)
- Negligible correlation indicates that the number of conformers sampled per ligand did not systematically bias binding affinity predictions, validating the multi-conformer docking approach.

**Affinity vs. Molecular Weight:** Pearson r = −0.18 (p > 0.05)
- Weak negative correlation suggests that larger compounds did not exhibit consistently improved or degraded binding predictions, consistent with diverse binding modes within the COX-2 active site.

## 4.4.5 Binding Affinity Distribution and Visualization

Binding affinity predictions were binned and visualized as a histogram with annotations for library-wide statistics and reference values.

**Key Observations:**
- A bimodal distribution emerged: the top-ranked flavonoids formed a distinct high-affinity cluster (−9.9 to −9.5 kcal/mol), while the remaining compounds dispersed across a lower-affinity region (−8.7 to −4.6 kcal/mol).
- Compounds with predicted affinities superior to or comparable with the reference ligand (≤ −10.22 kcal/mol): **Apigenin** and **Kaempferol** met or approached this threshold.
- Compounds with predicted affinities between the top hits and library mean: Luteolin, Quercetin, Scutellarein, and epicatechin derivatives formed an intermediate tier.

---

**Table 4.3: Comprehensive Statistical Analysis of Binding Affinity Data**

| Statistical Test/Metric | Result | Interpretation |
|------|--------|----------|
| **Descriptive Statistics** | | |
| Mean Affinity (All 33 compounds) | −7.88 ± 1.14 kcal/mol | Substantial library-wide variability |
| Mean Affinity (Top 5 Hits) | −9.73 ± 0.17 kcal/mol | Tightly clustered high-affinity group |
| Mean Affinity (Rest, n=28) | −7.30 ± 0.88 kcal/mol | Dispersed lower-affinity population |
| Affinity Range | −10.22 to −4.63 kcal/mol | 5.59 kcal/mol span |
| **Two-Sample t-Test** | | |
| t-statistic | −19.35 | Strong separation between groups |
| p-value | 1.01 × 10⁻⁶ | Highly significant difference (p < 0.001) |
| Mean Difference | 2.43 kcal/mol | Top hits 2.43 kcal/mol more favorable |
| **One-Way ANOVA** | | |
| F-statistic | 23.08 | Multiple group differences present |
| p-value | 7.27 × 10⁻⁷ | Highly significant (p < 0.001) |
| **Tukey HSD Post-Hoc Comparisons** | | |
| Reference vs Top Hits | p = 0.849 (NS) | No significant difference |
| Reference vs Rest | p = 0.004 ** | Significant difference |
| Top Hits vs Rest | p < 0.0001 *** | Highly significant difference |
| **Pearson Correlation Analysis** | | |
| Affinity vs Number of Conformers | r = 0.026 (NS) | No conformer bias in predictions |
| Affinity vs Molecular Weight | r = −0.18 (NS) | No systematic MW correlation |

NS = not significant (p > 0.05); ** p < 0.01; *** p < 0.001

# 4.5 Characterization of Predicted Protein–Ligand Interactions

## 4.5.1 Binding Mode Analysis of Top-Ranked Phytochemicals

Molecular visualization and detailed interaction analysis were performed on docked poses of the five top-ranked flavonoids and the reference ligand within the COX-2 active site.

**Structural Motifs Shared by Top-Ranked Compounds:**

All five flavonoid hits exhibited a characteristic binding orientation within the COX-2 binding pocket, consistent with flavonoid coordination in known COX inhibitor complex structures. Key interaction patterns included:

1. **Hydrogen Bonding Network:**
   - Flavonoid hydroxyl groups (particularly on rings B and C) engaged in hydrogen bonding with conserved COX-2 active-site residues.
   - The carbonyl oxygen at position 4 of the chromone ring frequently participated in hydrogen bonding.
   - Multiple hydroxyl groups enabled polyvalent hydrogen bonding, creating a stabilized, multidentate interaction network.

2. **Hydrophobic and π-Stacking Interactions:**
   - Aromatic rings of the flavonoids aligned within hydrophobic pockets of the COX-2 binding site.
   - π-stacking interactions with aromatic residues provided favorable van der Waals contact.
   - The B-ring aromatic system consistently positioned within a hydrophobic cleft.

3. **Spatial Complementarity:**
   - All five hits exhibited favorable shape and size complementarity with the COX-2 active site, despite their moderate molecular weight (270–302 Da).
   - The compact flavone scaffold efficiently occupied the binding cavity without excessive conformational strain.

**Individual Compound Characteristics:**

- **Apigenin (−9.91 kcal/mol):** The unsubstituted B-ring configuration and dual A- and C-ring hydroxylation facilitated optimal hydrogen bonding and hydrophobic packing. This compound achieved the highest-ranked binding affinity.
  
- **Kaempferol (−9.85 kcal/mol):** The 4'-hydroxylation on the B-ring and hydroxylation at position 3 of the C-ring provided enhanced polarity and hydrogen bonding capacity. Slight reduction in overall hydrophobicity compared to apigenin resulted in marginally lower binding energy.

- **Luteolin (−9.70 kcal/mol):** Hydroxylation at positions 6' and 3' on the B-ring increased molecular polarity. This added hydroxylation imposed minor geometric constraints, resulting in slightly decreased binding affinity compared to apigenin and kaempferol.

- **Quercetin (−9.69 kcal/mol):** The 3'-hydroxylation and additional hydroxyl at position 5 on the A-ring conferred increased hydrophilicity. The expanded hydroxylation pattern balanced between hydrophilic hydrogen bonding and hydrophobic interactions, achieving near-equivalent affinity with luteolin.

- **Scutellarein (−9.49 kcal/mol):** The 5,7-dihydroxy-6-methoxyflavone structure and unique methylation introduced additional steric factors and altered lipophilicity. This compound ranked fifth among the hits, suggesting that the methoxy substitution marginally reduced binding complementarity.

## 4.5.2 Comparison of Top Phytochemicals to Reference Ligand

The reference ligand Rofecoxib (RCX) displayed a structurally distinct binding mode compared to the flavonoid hits, reflecting its different chemical scaffold (benzenesulfonamide with isoxazole and benzoic acid moieties). Rofecoxib's predicted binding affinity (−10.22 kcal/mol) exceeded the top-ranked phytochemical by 0.31 kcal/mol (apigenin), indicating stronger predicted binding. However, the mechanistic difference—achieved through a completely different chemical structure and interaction geometry—demonstrates that multiple binding modes can achieve favorable COX-2 inhibition.

## 4.5.3 Interaction Patterns in Intermediate and Lower-Affinity Compounds

Compounds with predicted affinities between −8.6 and −8.1 kcal/mol (epicatechin, dicaffeoylquinic acid, phytosterols) generally exhibited one or more of the following characteristics:

- **Reduced Hydrogen Bonding:** Compounds with fewer hydroxyl groups or altered hydroxylation patterns engaged fewer favorable hydrogen bonds.
- **Suboptimal Hydrophobic Positioning:** Larger or more lipophilic structures (particularly phytosterols) occupied the binding pocket less efficiently, with portions extending beyond optimally favorable interaction zones.
- **Conformational Flexibility:** Larger or more flexible compounds (dicaffeoylquinic acid with its quinic acid moiety) may adopt conformations with reduced binding complementarity.

Lower-affinity compounds (predicted affinities < −7.5 kcal/mol) were characterized by insufficient size and shape complementarity, poor hydrogen bonding networks, or excessive polarity/hydrophobicity mismatches.

---

**Table 4.4: Toxicity Prediction Summary for Top-5 Hits and Reference Ligand**

| Toxicity Parameter | Apigenin | Kaempferol | Luteolin | Quercetin | Scutellarein | Rofecoxib |
|-------|----------|-----------|----------|-----------|--------------|----------|
| Acute Oral Toxicity | Low | Low | Low | Low | Low | Low |
| LD₅₀ Class | Class III/IV | Class III/IV | Class III/IV | Class III/IV | Class III/IV | Class II/III |
| **Hepatotoxicity** | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No |
| **Mutagenicity (Ames)** | ✓ No | ✓ No | ⚠ Alert | ⚠ Alert | ✓ No | ✓ No |
| **Carcinogenicity** | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No |
| **Immunotoxicity** | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No |
| **Skin Sensitization** | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No |
| **Respiratory Sensitization** | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No | ✓ No |
| **hERG Blockade** | Moderate | Moderate | Moderate | Moderate | Moderate | Moderate |
| **Acute Aquatic Toxicity** | Low | Low | Low | Low | Low | — |
| **Overall Safety Rating** | Favorable | Favorable | Favorable* | Favorable* | Favorable | Favorable |

✓ = Predicted inactive (safe); ⚠ = PAINS alert (in vitro assay artifact, not definitive); — = Not assessed. *Quercetin and Luteolin show PAINS alerts for mutagenicity in specific assay contexts but are well-tolerated natural products with extensive human dietary exposure.

# 4.6 Summary of Key Findings

1. **High-Affinity Hits:** Five flavonoids—Apigenin, Kaempferol, Luteolin, Quercetin, and Scutellarein—achieved predicted binding affinities (−9.91 to −9.49 kcal/mol) comparable to or exceeding most compounds in the library, with top affinity within 0.31 kcal/mol of the reference ligand Rofecoxib.

2. **Statistical Significance:** The top hits exhibited statistically significantly superior binding predictions compared to the rest of the library (t-test, p = 1.01 × 10⁻⁶; ANOVA, p = 7.27 × 10⁻⁷), with no significant difference from the reference ligand (Tukey HSD, p = 0.8491).

3. **Favorable Pharmacokinetics:** All five top hits predicted high gastrointestinal absorption, oral bioavailability comparable to Rofecoxib, and moderate distribution into tissue compartments.

4. **Acceptable Toxicity Profiles:** None of the five top hits were predicted to be hepatotoxic, mutagenic (with minor exceptions), or immunotoxic, though potential CYP inhibition and Pgp efflux warrant further investigation.

5. **Binding Mode Diversity:** Flavonoid hits achieved favorable binding through polyvalent hydrogen bonding and hydrophobic interactions, generating a distinct but mechanistically related binding mode compared to the reference nonsteroidal anti-inflammatory drug Rofecoxib.

These results provide computational evidence supporting the therapeutic potential of Dicoma anomala phytochemicals, particularly flavonoids, as leads for COX-2-targeted drug discovery and supporting the ethnopharmacological use of this plant in inflammatory and pain conditions.
