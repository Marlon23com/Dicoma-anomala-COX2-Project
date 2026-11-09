# Dicoma anomala COX-2 Screening Project

This repository contains a computational drug-discovery workflow evaluating phytochemicals from Dicoma anomala as potential Cyclooxygenase-2 (COX-2) inhibitors for anti-inflammatory applications.

## Project overview

The project combines:

- literature mining and manual curation of Dicoma anomala compounds
- ligand preparation and docking workflows
- COX-2 protein preparation and receptor docking
- ADMET and toxicity screening
- statistical analysis of docking outcomes
- visualization of molecular interactions and key findings

The study focuses on identifying natural-product hits with favorable binding affinity and acceptable pharmacokinetic profiles when compared with the reference ligand Rofecoxib (RCX).

## Repository structure

- `00_Documentation/` – project notes, literature summaries, proposal and results chapters
- `01_Literature_Data/` – curated compound data and mining outputs
- `02_Ligands/` – ligand preparation and source files
- `03_Protein/` – protein preparation and receptor structures
- `04_ADMET/` – SwissADME and ProTox analyses
- `05_Docking/` – docking scripts, outputs and summary statistics
- `06_Visualization/` – molecular interaction and publication images
- `07_Analysis/` – downstream analysis and statistical assets
- `AutoDock-Vina/` – local AutoDock Vina source and documentation
- `scripts/` – supporting scripts and environment files

## Key findings

The screening workflow identified several flavonoid compounds, including Apigenin, Kaempferol, Luteolin, Quercetin, and Scutellarein, as top-ranked COX-2 binders with predicted affinities close to the reference inhibitor. These compounds were further assessed for ADMET properties, binding behavior, and statistical significance.

## Notes

This repository intentionally excludes the thesis draft materials in `00_Documentation/thesis_drafts/`.

## Reproducibility

The project includes environment and dependency files to support reproduction of the computational workflow.
