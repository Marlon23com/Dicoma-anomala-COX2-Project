# TODO: This code to be reviewed by Claude Code for validation and optimization

import os
import sys
from pathlib import Path
import logging

# 1. Force the use of openbabel-wheel's internal data
# Find the location of the openbabel module
try:
    import openbabel
    ob_path = Path(openbabel.__file__).parent
    # The data folder is usually inside the openbabel package directory in the wheel
    data_dir = ob_path / "data"
    
    if data_dir.exists():
        os.environ['BABEL_DATADIR'] = str(data_dir)
        print(f"Using Open Babel data from: {data_dir}")
    else:
        print("Warning: Could not find internal data dir. Trying system default.")
except ImportError:
    print("Error: openbabel module not found. Ensure 'openbabel-wheel' is installed in this environment.")
    sys.exit(1)

# Now import pybel after setting the env var
from openbabel import pybel

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def prepare_protein_fixed(input_pdb: str, output_pdbqt: str):
    """
    Prepares COX-2 protein for docking using openbabel-wheel with fixed paths.
    """
    input_path = Path(input_pdb)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_pdb}")
        return False

    try:
        logger.info("Loading PDB file...")
        mol = next(pybel.readfile("pdb", str(input_path)))
        
        logger.info("Removing water and heteroatoms...")
        keep_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        atoms_to_delete = []
        for atom in mol.atoms:
            res_name = atom.OBAtom.GetResidue().GetName().strip()
            if res_name not in keep_residues:
                atoms_to_delete.append(atom.idx)
                
        for idx in sorted(atoms_to_delete, reverse=True):
            mol.OBMol.DeleteAtom(mol.OBMol.GetAtom(idx))
            
        logger.info("Adding polar hydrogens...")
        mol.addh()
        
        logger.info("Assigning Gasteiger charges...")
        mol.calccharges()
        
        logger.info("Saving as PDBQT...")
        mol.write("pdbqt", str(output_pdbqt), overwrite=True)
        
        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size > 1000:
            size_kb = Path(output_pdbqt).stat().st_size / 1024
            logger.info(f"✅ Successfully created: {output_pdbqt} ({size_kb:.1f} KB)")
            return True
        else:
            logger.error("❌ Output file is missing or too small.")
            return False
            
    except Exception as e:
        logger.error(f"❌ Error during preparation: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    INPUT_PDB = "./raw_pdb/5KIR.pdb"
    OUTPUT_PDBQT = "./prepared_pdbqt/5KIR_prepared.pdbqt"
    
    Path("./prepared_pdbqt").mkdir(parents=True, exist_ok=True)
    
    success = prepare_protein_fixed(INPUT_PDB, OUTPUT_PDBQT)
    
    if success:
        print("\n🎉 Protein preparation complete! You are ready for docking.")
    else:
        print("\n💥 Protein preparation failed. Check logs above.")