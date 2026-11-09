# TODO: This code to be reviewed by Claude Code for validation and optimization

import subprocess
from pathlib import Path
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def prepare_protein_direct(input_pdb: str, output_pdbqt: str):
    """
    Converts PDB directly to PDBQT using Open Babel 2.4.1.
    Handles deletion of water/ligand, addition of hydrogens, and charges in one command.
    """
    input_path = Path(input_pdb)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_pdb}")
        return False

    # Define the exact path to your Open Babel executable
    obabel_exe = r"C:\Program Files\OpenBabel-2.4.1\obabel.exe"
    
    if not os.path.exists(obabel_exe):
        logger.error(f"Open Babel executable not found at: {obabel_exe}")
        return False

    logger.info("Starting direct conversion: PDB -> PDBQT...")
    logger.info("Actions: Delete Water/Ligand, Add H, Assign Gasteiger Charges.")
    
    # Command breakdown:
    # input_pdb : The source file
    # -O output_pdbqt : The destination file
    # --delete HOH : Remove water molecules
    # --delete RCX : Remove Rofecoxib (the ligand)
    # -h : Add polar hydrogens appropriate for pH 7.4
    # --partialcharge gasteiger : Assign Gasteiger charges required for Vina
    cmd = [
        obabel_exe, 
        str(input_path), 
        "-O", output_pdbqt,
        "--delete", "HOH",
        "--delete", "RCX",
        "-h",
        "--partialcharge", "gasteiger"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if Path(output_pdbqt).exists():
            size_kb = Path(output_pdbqt).stat().st_size / 1024
            if size_kb > 100: 
                logger.info(f" Successfully created: {output_pdbqt} ({size_kb:.1f} KB)")
                return True
            else:
                logger.error(f" Conversion failed: Output PDBQT is too small ({size_kb:.1f} KB).")
                if result.stderr:
                    logger.error(f"Open Babel Error Log:\n{result.stderr}")
                return False
        else:
            logger.error(" Conversion failed: Output PDBQT file was not created.")
            if result.stderr:
                logger.error(f"Open Babel Error Log:\n{result.stderr}")
            return False
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Open Babel execution error: {e.stderr}")
        return False

if __name__ == "__main__":
    # Paths
    INPUT_PDB = "./raw_pdb/5KIR.pdb"
    OUTPUT_PDBQT = "./prepared_pdbqt/5KIR_prepared.pdbqt"
    
    # Ensure output directory exists
    Path("./prepared_pdbqt").mkdir(parents=True, exist_ok=True)
    
    success = prepare_protein_direct(INPUT_PDB, OUTPUT_PDBQT)
    
    if success:
        print("\n Protein preparation complete! You are ready for docking.")
    else:
        print("\n Protein preparation failed. Check logs above.")


