# TODO: This code to be reviewed by Claude Code for validation and optimization

import os
import subprocess
import glob
import pandas as pd
from pathlib import Path
import logging
import time

# Setup Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_vina_docking(receptor: str, ligand: str, output: str, config: dict):
    """
    Runs AutoDock Vina for a single ligand conformer.
    Returns the best binding affinity (kcal/mol).
    """
    cmd = [
        "vina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--out", output,
        "--center_x", str(config['center_x']),
        "--center_y", str(config['center_y']),
        "--center_z", str(config['center_z']),
        "--size_x", str(config['size_x']),
        "--size_y", str(config['size_y']),
        "--size_z", str(config['size_z']),
        "--exhaustiveness", str(config['exhaustiveness']),
        "--num_modes", str(config['num_modes']),
        "--energy_range", str(config['energy_range'])
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse affinity from stdout
        # Vina output usually looks like: "   1     -8.5      0.000      0.000      0.000"
        for line in result.stdout.splitlines():
            if "Affinity" in line or (line.strip() and line.strip()[0].isdigit()):
                parts = line.split()
                # Find the first float that looks like an energy (negative number)
                for part in parts:
                    try:
                        val = float(part)
                        if val < 0: 
                            return val
                    except ValueError:
                        continue
        return None
    except subprocess.CalledProcessError as e:
        logger.error(f"Vina failed for {os.path.basename(ligand)}: {e.stderr}")
        return None

def main():
    # --- CONFIGURATION ---
    RECEPTOR = "03_Protein/prepared_pdbqt/5KIR_prepared_clean.pdbqt"
    LIGAND_DIR = "02_Ligands/raw_sdf/PDBQT_20Confs"
    OUTPUT_DIR = "05_Docking/results_multi_conf"
    
    # Grid Box Parameters (Centered on RCX)
    CONFIG = {
        'center_x': 22.5,
        'center_y': 1.0,
        'center_z': 36.5,
        'size_x': 25,
        'size_y': 25,
        'size_z': 25,
        'exhaustiveness': 32,
        'num_modes': 10,
        'energy_range': 4.0
    }
    
    # Create output directory
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    
    # Get all PDBQT files
    all_files = glob.glob(os.path.join(LIGAND_DIR, "*.pdbqt"))
    
    if not all_files:
        logger.error(f"No PDBQT files found in {LIGAND_DIR}")
        return

    logger.info(f"Found {len(all_files)} PDBQT files to dock.")
    
    results = []
    start_time = time.time()
    
    for i, ligand_file in enumerate(all_files):
        ligand_name = os.path.basename(ligand_file)
        logger.info(f"[{i+1}/{len(all_files)}] Docking: {ligand_name}")
        
        # Define output path
        out_file = os.path.join(OUTPUT_DIR, ligand_name.replace('.pdbqt', '_out.pdbqt'))
        
        # Run Docking
        affinity = run_vina_docking(RECEPTOR, ligand_file, out_file, CONFIG)
        
        if affinity is not None:
            logger.info(f"  -> Affinity: {affinity} kcal/mol")
            results.append({
                'Filename': ligand_name,
                'Affinity_kcal_mol': affinity,
                'Output_File': out_file
            })
        else:
            logger.warning(f"  -> Failed to get affinity for {ligand_name}")
            
    # Save Results
    if results:
        df = pd.DataFrame(results)
        # Sort by best affinity (most negative first)
        df = df.sort_values(by='Affinity_kcal_mol', ascending=True)
        df.to_csv("05_Docking/docking_results_summary.csv", index=False)
        logger.info(f"\n✅ Docking Complete! Results saved to 05_Docking/docking_results_summary.csv")
        
        # Print Top 5
        print("\n--- Top 5 Best Binders ---")
        print(df.head(5).to_string(index=False))
    else:
        logger.error("No successful docking results to save.")

    end_time = time.time()
    logger.info(f"Total Time: {(end_time - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()