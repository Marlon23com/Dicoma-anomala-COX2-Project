# TODO: This code to be reviewed by Claude Code for validation and optimization

import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def evaluate_admet_hits(input_file, output_file, reference_compound="Rofecoxib"):
    """
    Evaluates ADMET properties based on Table 3.1 thresholds and calculates 
    average normalized deviation to select hits for MD simulations.
    """
    try:
        # Load Data
        df = pd.read_csv(input_file)
        logger.info(f"Loaded {len(df)} compounds from {input_file}")
        
        # --- 1. Define Column Mapping ---
        # IMPORTANT: Update these keys to match the EXACT headers in your ADMETLab CSV
        # Check your CSV file to see if it uses 'Caco-2', 'caco2', 'Caco2', etc.
        col_map = {
            'Caco2': 'caco2',          # Absorption
            'HIA': 'hia',               # Human Intestinal Absorption (Probability)
            'PPB': 'PPB',               # Plasma Protein Binding (%)
            'Fu': 'Fu',                 # Fraction Unbound (%)
            'CYP3A4_sub': 'CYP3A4-sub', # Substrate Probability
            'CYP3A4_inh': 'CYP3A4-inh', # Inhibitor Probability
            'Cl_plasma': 'cl-plasma',   # Clearance
            't_half': 't0.5',           # Half-life
            'Ames': 'Ames',             # Mutagenicity Probability
            'hERG': 'hERG'              # Cardiotoxicity Probability
        }

        # --- 2. Define Thresholds (From Table 3.1) ---
        # Direction: 'high_is_good' means Value > Threshold is positive deviation.
        # Direction: 'low_is_good' means Value < Threshold is positive deviation.
        criteria = {
            'Caco2': {'threshold': -5.15, 'direction': 'high_is_good'}, # >= -5.15
            'HIA': {'threshold': 0.7, 'direction': 'high_is_good'},     # >= 0.7 (High Absorption is good)
            'PPB': {'threshold': 90, 'direction': 'low_is_good'},       # <= 90%
            'Fu': {'threshold': 5, 'direction': 'high_is_good'},        # >= 5%
            'CYP3A4_sub': {'threshold': 0.3, 'direction': 'low_is_good'}, # <= 0.3
            'CYP3A4_inh': {'threshold': 0.3, 'direction': 'low_is_good'}, # <= 0.3
            'Cl_plasma': {'threshold': 5, 'direction': 'low_is_good'},  # < 5
            't_half': {'threshold': 0.7, 'direction': 'low_is_good'},   # <= 0.7
            'Ames': {'threshold': 0.7, 'direction': 'low_is_good'},     # <= 0.7 (Non-toxic)
            'hERG': {'threshold': 0.7, 'direction': 'low_is_good'}      # <= 0.7 (Non-toxic)
        }

        # --- 3. Calculate Deviations ---
        deviations = []
        
        for index, row in df.iterrows():
            # Try to find the compound name column
            compound_name = row.get('Compound Name', row.get('name', row.get('Molecule', f'Index_{index}')))
            compound_devs = {}
            valid_count = 0
            
            for prop, config in criteria.items():
                col_name = col_map.get(prop)
                if col_name not in df.columns:
                    # Try case-insensitive search if exact match fails
                    matching_cols = [c for c in df.columns if c.lower() == col_name.lower()]
                    if matching_cols:
                        col_name = matching_cols[0]
                    else:
                        continue
                
                val = row[col_name]
                try:
                    val = float(val)
                except:
                    continue
                
                threshold = config['threshold']
                
                # Calculate Deviation
                if config['direction'] == 'high_is_good':
                    dev = val - threshold
                else:
                    dev = threshold - val
                    
                compound_devs[prop] = dev
                valid_count += 1
                
            # Calculate Average Normalized Deviation
            if valid_count > 0:
                avg_dev = np.mean(list(compound_devs.values()))
            else:
                avg_dev = 0
                
            deviations.append({
                'Compound': compound_name,
                'Avg_Deviation': avg_dev,
                **compound_devs
            })
            
        dev_df = pd.DataFrame(deviations)
        
        # --- 4. Select Hits Based on Reference (Rofecoxib) ---
        ref_score_row = dev_df[dev_df['Compound'] == reference_compound]
        
        if not ref_score_row.empty:
            ref_score = ref_score_row['Avg_Deviation'].values[0]
            logger.info(f"Reference Compound ({reference_compound}) Avg Deviation: {ref_score:.4f}")
            
            # Select compounds with HIGHER deviation (better profile) than reference
            selected_hits = dev_df[dev_df['Avg_Deviation'] > ref_score]
            
            if selected_hits.empty:
                logger.info("No compounds performed better than Rofecoxib. Selecting top 3 closest.")
                # Sort by deviation descending and take top 3 excluding reference
                non_ref = dev_df[dev_df['Compound'] != reference_compound]
                selected_hits = non_ref.nlargest(3, 'Avg_Deviation')
        else:
            logger.warning(f"Reference compound '{reference_compound}' not found. Selecting top 5 by deviation.")
            selected_hits = dev_df.nlargest(5, 'Avg_Deviation')
            
        # Sort selected hits
        selected_hits = selected_hits.sort_values(by='Avg_Deviation', ascending=False)
        
        # Save Results
        selected_hits.to_csv(output_file, index=False)
        logger.info(f"Selected {len(selected_hits)} hits for MD simulations. Saved to {output_file}")
        
        print("\n--- Selected Hits for MD Simulations ---")
        print(selected_hits[['Compound', 'Avg_Deviation']].to_string(index=False))
        
        return selected_hits

    except Exception as e:
        logger.error(f"Error in evaluation: {e}")
        raise

if __name__ == "__main__":
    # INPUT: Your ADMETLab CSV
    INPUT_CSV = "admetlab3_analysis_for_top5_and_reference.csv" 
    OUTPUT_CSV = "selected_hits_for_MD.csv"
    REFERENCE = "Rofecoxib" 
    
    evaluate_admet_hits(INPUT_CSV, OUTPUT_CSV, reference_compound=REFERENCE)