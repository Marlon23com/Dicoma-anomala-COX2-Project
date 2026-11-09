# CHAPTER 3: MATERIALS AND METHODS

This chapter describes the materials, computational tools, and analytical workflows used in the study. It outlines the procedures employed to identify, curate, prepare, and evaluate Dicoma anomala phytochemical ligands, to select and process the COX-2 target protein, to screen compounds for drug-likeness and ADMET properties, and to perform docking, interaction visualization, and cheminformatics analysis.

# 3.1 Research Design and Computational Workflow

This study adopted a materials and methods framework appropriate for a fully in silico investigation of Dicoma anomala phytochemicals as potential Cyclooxygenase-2 (COX-2) inhibitors. The computational workflow was structured as a sequence of interdependent phases that moved from compound identification to ligand preparation, protein selection, ADMET screening, molecular docking, interaction visualization, cheminformatics profiling, and statistical analysis. The methodology emphasised reproducibility through the use of documented software versions and standardised file formats, with scripting performed in Python 3.12.3 and data processing conducted using scientific libraries such as pandas, numpy, scipy, statsmodels, matplotlib, and seaborn.

# 3.2 Identification and Preparation of Phytochemical Ligands

## 3.2.1 Phytochemical Compilation and Curation

A targeted library of phytochemicals associated with Dicoma anomala was assembled using a hybrid combination of automated literature mining and manual curation. Custom Python scripts were developed to query scientific literature sources and chemical databases, including Scopus, Web of Science, PubMed, and the PubChem Power User Gateway (PUG-REST). Search strategies employed controlled vocabularies and Medical Subject Heading (MeSH) terms such as "Dicoma anomala," "phytochemical composition," "secondary metabolites," and "anti-inflammatory compounds." Query results were retrieved in structured form, and error handling was implemented to manage API rate limits, malformed responses, and inconsistent metadata.

Returned records were programmatically filtered to retain only compounds explicitly linked to Dicoma anomala at the species level, with taxonomic validation performed through comparison against authoritative botanical references and phytochemical review literature. Duplicate entries were identified and collapsed using canonical identifiers, including SMILES strings and InChIKeys, to produce a non-redundant compound list. Manual verification followed automated retrieval to confirm botanical origin, excluding ambiguous entries attributed only to related Dicoma species unless corroborated by cross-referenced phytochemical reviews. Selected compounds were exported in Structure Data File (SDF) format to preserve stereochemical information, and canonical SMILES strings were extracted in parallel to support subsequent cheminformatics processing.

## 3.2.2 Ligand Preparation

Ligand structures were retrieved in SMILES format from PubChem and then converted to three-dimensional representations for docking preparation. Open Babel (version 3.1.1) was used for initial format conversion and to produce 3D SDF files from SMILES, preserving stereochemical configuration where present. Additional conformer generation and geometry optimisation were performed using cheminformatics routines, including the ETKDGv3 algorithm for 3D coordinate sampling and the MMFF94 force field for energy minimisation. This sequence of steps reduced steric strain and produced plausible low-energy conformers suitable for docking.

Prepared ligand structures were further processed by adding polar hydrogen atoms at physiological pH (7.4) and assigning Gasteiger partial charges, ensuring compatibility with AutoDock Vina's scoring function. The final ligand files were converted to PDBQT format using Open Babel, retaining atom types, charges, and torsional degrees of freedom required for docking.

A parallelised Python script was employed to automate ligand preparation across the compound library. The script handled file conversion, conformer sampling, energy minimisation, charge assignment, and PDBQT generation in batch mode, thereby improving computational efficiency and reproducibility.

# 3.3 Target Protein Selection and Preparation

## 3.3.1 Protein Retrieval

The target protein for molecular docking was human Cyclooxygenase-2 (COX-2) in complex with the inhibitor Rofecoxib, obtained from the Protein Data Bank (PDB ID: 5KIR). This structure was chosen for its available co-crystallised ligand and well-defined binding site.

## 3.3.2 Protein Preparation

Protein preparation was carried out using UCSF Chimera (version 1.11.2) and AutoDock Tools. The PDB entry was preprocessed to remove non-essential heteroatoms, including crystallographic water molecules and ions not directly involved in ligand binding, while retaining active-site residues and protein backbone integrity.

Polar hydrogens were added to the protein to reflect physiological protonation states, and Gasteiger-Marsili partial charges were assigned to all atoms within the binding site. The prepared protein structure was saved in PDBQT format to encode atomic coordinates, partial charges, and AutoDock-compatible atom types. The coordinate frame was preserved for grid definition and docking execution.

# 3.4 In Silico ADMET and Drug-Likeness Profiling

## 3.4.1 Drug-Likeness Assessment

Drug-likeness evaluation was performed using SwissADME, a web-based tool for physicochemical and lipophilicity analysis. For each compound, the following descriptors were calculated: molecular weight, octanol-water partition coefficient (logP), topological polar surface area (TPSA), hydrogen bond donors, hydrogen bond acceptors, and rotatable bond count. These descriptors were assessed against established drug-likeness rules, including Lipinski's Rule of Five, to identify compounds with properties indicative of acceptable oral bioavailability.

## 3.4.2 ADMET Prediction

ADMET prediction involved in silico evaluation of absorption, distribution, metabolism, excretion, and toxicity properties. Publicly available prediction platforms such as ProTox-II were used to estimate acute oral toxicity, hepatotoxicity, immunotoxicity, and mutagenicity. Compounds with predicted liabilities in critical ADMET domains were noted for downstream interpretation and prioritisation during interaction analysis.

# 3.5 Molecular Docking Protocol

## 3.5.1 Docking Software and Algorithm

Molecular docking was conducted using AutoDock Vina (version 1.2.3), which employs a stochastic global search coupled with an empirical scoring function. Docking simulations treated the protein as rigid and allowed ligand flexibility through rotatable bonds defined in the PDBQT files.

## 3.5.2 Grid Box Parameters

The docking grid box was defined around the COX-2 active site to encompass the co-crystallised ligand-binding pocket. The box dimensions were set to 20 Å × 20 Å × 20 Å, centred on the coordinates of the Rofecoxib ligand in PDB 5KIR, to ensure the search space captured the known inhibitor binding location and accommodated the size of the prepared phytochemicals.

## 3.5.3 Docking Execution

Docking execution was automated using custom Python scripts, which managed the sequential launching of AutoDock Vina jobs and the collection of output files. Parameters included exhaustiveness = 8 to control the depth of conformational sampling, num_modes = 10 to limit the maximum number of output poses per ligand, and energy_range = 3.0 kcal/mol to define the energy clustering window. Batch processing was performed across available CPU cores to increase computational efficiency.

All docking outputs were generated in PDBQT format, and file naming conventions encoded both ligand identifiers and conformer indices to preserve traceability.

# 3.6 Interaction Analysis and Visualization

## 3.6.1 Visualization Tools

Interaction analysis and visualization were performed using molecular graphics tools such as UCSF Chimera and PyMOL. These tools were used to inspect docked ligand poses, visualise protein-ligand complexes, and create publication-quality figures illustrating binding modes.

## 3.6.2 Interaction Criteria

Protein-ligand interactions were evaluated based on standard criteria, including hydrogen bonding, hydrophobic contacts, pi stacking, and electrostatic interactions. Distance cutoffs and geometric criteria recommended in structural biology literature were applied to identify relevant interactions within the COX-2 active site.

# 3.7 Statistical and Cheminformatics Analysis

## 3.7.1 Cheminformatics Profiling

Cheminformatics profiling was conducted using Python scripting and relevant libraries to calculate molecular descriptors, generate canonical SMILES, and compute conformer counts. Where available, RDKit was employed to generate three-dimensional conformers, compute molecular fingerprints, and standardise chemical representations. Data handling and descriptor aggregation were performed using pandas and numpy.

## 3.7.2 Statistical Analysis

Statistical analysis was performed on docking-derived affinity predictions and descriptor data using scipy and statsmodels. The workflow included group comparisons, distribution analysis, and correlation testing to investigate relationships between ligand properties and predicted binding affinity. Visualization of statistical distributions was performed with matplotlib and seaborn to support the interpretation of the computational results in later chapters.

# 3.8 Ethical Considerations

This research was purely computational and did not involve human participants, animals, or primary biological samples. All data were obtained from publicly accessible databases and published literature. Ethical considerations focused on compliance with database usage terms, attribution of source materials, and transparent reporting of methods. No ethical approval was required for the computational analyses described in this chapter.

