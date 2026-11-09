# TODO: This code to be reviewed by Claude Code for validation and optimization

import os
import subprocess
import logging
import sys
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops

# --- Configuration ---
input_dir = "."  # The directory where the script is located (raw_sdf)
output_dir = os.path.join(input_dir, "PDBQT_20Confs")  # Output directory
log_file = "conversion_20confs.log"
max_parallel_jobs = 4  # Reduced for stability on Windows
num_conformers = 20    # Number of conformers per molecule
charge_model = "gasteiger"

# HARDCODED PATH TO OBABEL (Change this if your path is different)
# Based on your previous logs, this is where Conda installed it:
OBABEL_EXE = r"C:\Users\hp\Anaconda3\Library\bin\obabel.exe"

# Set up logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s', handlers=[
    logging.FileHandler(log_file),
    logging.StreamHandler()
])

def validate_obabel():
    if not os.path.exists(OBABEL_EXE):
        logging.error(f"Open Babel not found at: {OBABEL_EXE}")
        sys.exit(1)
    logging.info(f"Using Open Babel at: {OBABEL_EXE}")

def validate_input_dir():
    if not os.path.isdir(input_dir):
        logging.error(f"Input directory '{input_dir}' does not exist")
        sys.exit(1)
    sdf_files = [f for f in os.listdir(input_dir) if f.endswith(".sdf")]
    if not sdf_files:
        logging.error(f"No SDF files found in '{input_dir}'")
        sys.exit(1)
    else:
        logging.info(f"Found {len(sdf_files)} SDF files in '{input_dir}'")

def validate_output_dir():
    if not os.path.isdir(output_dir):
        logging.info(f"Creating output directory '{output_dir}'")
        os.makedirs(output_dir, exist_ok=True)

def process_molecule(args):
    sdf_file, base_name, idx, mol = args
    
    if mol is None:
        logging.error(f"Skipping invalid molecule in '{base_name}' at index {idx}")
        return

    try:
        # Add hydrogens and perceive aromaticity
        mol = Chem.AddHs(mol)
        rdmolops.Kekulize(mol)

        # Generate multiple conformers
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, pruneRmsThresh=0.5, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        if not cids:
            logging.error(f"Error generating conformers for molecule '{base_name}' at index {idx}")
            return

        for cid in cids:
            # Minimize energy using MMFF94
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)

            # Write to temporary PDB file
            temp_pdb = os.path.join(output_dir, f"{base_name}_{idx}_{cid}.pdb")
            Chem.MolToPDBFile(mol, temp_pdb, confId=cid)

            # Convert PDB to PDBQT using Open Babel
            pdbqt_file = os.path.join(output_dir, f"{base_name}_{idx}_{cid}.pdbqt")
            
            # Use list format for subprocess to avoid Windows shell parsing issues
            cmd = [
                OBABEL_EXE, 
                temp_pdb, 
                "-O", pdbqt_file, 
                "--partialcharge", charge_model
            ]
            
            # Run subprocess
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            if os.path.exists(temp_pdb):
                os.remove(temp_pdb)  # Remove temporary PDB file

            if result.returncode == 0:
                if os.path.exists(pdbqt_file) and os.path.getsize(pdbqt_file) > 0:
                    logging.info(f"SUCCESS: Converted conformer {cid} of {base_name}")
                else:
                    logging.error(f"EMPTY: Output file '{pdbqt_file}' is empty")
            else:
                err_msg = result.stderr.decode('utf-8', errors='ignore')
                logging.error(f"ERROR: Failed to convert conformer {cid} of {base_name}: {err_msg}")

    except RuntimeError as e:
        logging.error(f"Runtime error processing molecule {idx} in '{base_name}': {e}")
    except Exception as e:
        logging.error(f"Error processing molecule {idx} in '{base_name}': {e}")

def prepare_args(sdf_file):
    base_name = os.path.basename(sdf_file).replace('.sdf', '')
    
    if not os.path.isfile(sdf_file) or os.path.getsize(sdf_file) == 0:
        return []

    try:
        suppl = Chem.SDMolSupplier(sdf_file)
        args_list = []
        for idx, mol in enumerate(suppl):
            args_list.append((sdf_file, base_name, idx, mol))
        return args_list
    except Exception as e:
        logging.error(f"Error reading '{sdf_file}': {e}")
        return []

def main():
    # Initialize log file
    with open(log_file, 'w') as f:
        f.write(f"=== Conversion Started ==={os.linesep}")

    # Perform validations
    validate_obabel()
    validate_input_dir()
    validate_output_dir()

    # Get all SDF files in the input directory
    sdf_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.sdf')]

    # Prepare all arguments
    all_args = []
    for sdf_file in sdf_files:
        all_args.extend(prepare_args(sdf_file))

    if not all_args:
        logging.error("No valid molecules found to process.")
        return

    # Process files in parallel
    logging.info(f"Starting parallel conversion with max {max_parallel_jobs} jobs for {len(all_args)} conformers")
    
    with Pool(min(max_parallel_jobs, cpu_count())) as pool:
        pool.map(process_molecule, all_args, chunksize=1)

    # Log summary
    logging.info(f"SUMMARY: Conversion completed. Check '{output_dir}' for results.")

if __name__ == "__main__":
    main()