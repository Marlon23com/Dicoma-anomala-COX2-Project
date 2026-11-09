# TODO: This code to be reviewed by Claude Code for validation and optimization

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def validate_and_clean_sdf(input_dir: str, output_dir: str):
    """
    Validates SDF files, adds hydrogens, minimizes energy, and saves clean PDBQT-ready SDFs.
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Look for .sdf files specifically
    sdf_files = list(input_path.glob("*.sdf"))
    
    # If no .sdf found, check if files exist without extension or with .txt
    if not sdf_files:
        all_files = list(input_path.iterdir())
        logger.warning(f"No .sdf files found in {input_dir}. Found these files: {[f.name for f in all_files]}")
        return

    logger.info(f"Found {len(sdf_files)} SDF files to validate.")
    
    valid_count = 0
    failed_count = 0
    failed_list = []

    for sdf_file in sdf_files:
        try:
            # 1. Load Molecule
            supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
            mol = next(supplier)
            
            if mol is None:
                raise ValueError("Could not parse molecule")
            
            # 2. Check for Valid Structure
            if not mol.GetNumAtoms():
                raise ValueError("Empty molecule")
            
            # 3. Add Hydrogens (Crucial for Docking)
            mol = Chem.AddHs(mol)
            
            # 4. Embed/Generate 3D Coordinates (if missing or bad)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42  # Reproducibility
            res = AllChem.EmbedMolecule(mol, params)
            
            if res != 0:
                logger.warning(f"{sdf_file.name}: Failed to embed 3D coordinates. Trying UFF optimization on existing coords.")
                # Fallback: Just optimize whatever coords exist
                AllChem.UFFOptimizeMolecule(mol)
            else:
                # 5. Energy Minimization (UFF Force Field)
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            
            # 6. Save Cleaned SDF
            output_file = output_path / sdf_file.name
            writer = Chem.SDWriter(str(output_file))
            writer.write(mol)
            writer.close()
            
            valid_count += 1
            logger.info(f"✅ Validated: {sdf_file.name}")
            
        except Exception as e:
            failed_count += 1
            failed_list.append(sdf_file.name)
            logger.error(f"❌ Failed {sdf_file.name}: {str(e)}")
            
    # Summary
    logger.info("\n--- Validation Summary ---")
    logger.info(f"Valid: {valid_count}")
    logger.info(f"Failed: {failed_count}")
    if failed_list:
        logger.warning(f"Failed files: {failed_list}")
        with open(output_path / "validation_failures.txt", "w") as f:
            for name in failed_list:
                f.write(f"{name}\n")

if __name__ == "__main__":
    # Configure paths
    # Since the script is INSIDE raw_sdf, we look in the current directory "."
    INPUT_DIR = "." 
    # We save validated files to a new folder 'validated_sdf' inside raw_sdf
    OUTPUT_DIR = "./validated_sdf"
    
    validate_and_clean_sdf(INPUT_DIR, OUTPUT_DIR)