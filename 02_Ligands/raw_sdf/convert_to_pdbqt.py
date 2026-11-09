# TODO: This code to be reviewed by Claude Code for validation and optimization

import os
import subprocess
from pathlib import Path
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def convert_sdf_to_pdbqt(input_dir: str, output_dir: str):
    """
    Converts all SDF files in input_dir to PDBQT files in output_dir using Open Babel.
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    sdf_files = list(input_path.glob("*.sdf"))
    logger.info(f"Found {len(sdf_files)} SDF files to convert.")
    
    success_count = 0
    failed_count = 0
    
    for sdf_file in sdf_files:
        try:
            # Define output filename
            pdbqt_name = sdf_file.stem + ".pdbqt"
            pdbqt_file = output_path / pdbqt_name
            
            # Open Babel command
            # -h: Add hydrogens (already done, but safe to ensure)
            # -p 7.4: Assign protonation state at pH 7.4
            # --partialcharge gasteiger: Assign Gasteiger charges (Required for Vina)
            cmd = [
                "obabel", 
                str(sdf_file), 
                "-O", str(pdbqt_file), 
                "-p", "7.4", 
                "--partialcharge", "gasteiger"
            ]
            
            # Run Open Babel
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0 and pdbqt_file.exists():
                logger.info(f"✅ Converted: {sdf_file.name} -> {pdbqt_name}")
                success_count += 1
            else:
                logger.error(f"❌ Failed to convert {sdf_file.name}: {result.stderr}")
                failed_count += 1
                
        except Exception as e:
            logger.error(f"❌ Error processing {sdf_file.name}: {str(e)}")
            failed_count += 1
            
    logger.info("\n--- Conversion Summary ---")
    logger.info(f"Successful: {success_count}")
    logger.info(f"Failed: {failed_count}")

if __name__ == "__main__":
    # Input is the folder where validate_structures.py saved the clean SDFs
    INPUT_DIR = "./validated_sdf" 
    OUTPUT_DIR = "./prepared_pdbqt"
    
    convert_sdf_to_pdbqt(INPUT_DIR, OUTPUT_DIR)