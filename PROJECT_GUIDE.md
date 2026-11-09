# Dicoma anomala COX-2 Screening Project Guide

This dissertation project repository is organized to support end-to-end in silico screening of Dicoma anomala phytochemicals against COX-2 for sinusitis treatment.

## Folder overview

- `00_Documentation/`: Project documentation, literature review, proposal drafts, and thesis drafts.
- `01_Literature_Data/`: Literature mining results, compound lists, and curated references.
- `02_Ligands/`: Ligand source files, prepared docking inputs, and ligand preparation scripts.
- `03_Protein/`: Protein source structures, prepared receptor files, and active site definition scripts.
- `04_ADMET/`: ADMET screening results, toxicity predictions, and related processing scripts.
- `05_Docking/`: Docking configuration, control docking data, results summaries, and batch docking runs.
- `06_Visualization/`: Visualization protocols, interaction figures, and PyMOL session files.
- `07_Analysis/`: Cheminformatics, machine learning, and statistical analysis assets.
- `08_Backups/`: Backup snapshots organized by daily, weekly, and final submission folders.

## How to use each folder

1. Start in `01_Literature_Data/` to collect and curate compound candidates.
2. Prepare ligand files in `02_Ligands/` and protein files in `03_Protein/`.
3. Run ADMET filtering in `04_ADMET/` before docking.
4. Configure and execute docking workflows in `05_Docking/`.
5. Generate interaction visuals in `06_Visualization/` and analyze results in `07_Analysis/`.
6. Save major checkpoints in `08_Backups/`.

## Environment and dependencies

- Use `scripts/requirements.txt` for pip-based installs:
  - `pip install -r scripts/requirements.txt`
- Use `scripts/environment.yml` for Conda environment creation:
  - `conda env create -f scripts/environment.yml`
  - `conda activate dicoma_anomala_cox2`
