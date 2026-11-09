# clean_receptor.py
import os

input_file = "5KIR_prepared.pdbqt"
output_file = "5KIR_prepared_clean.pdbqt"

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        # Only write lines that are NOT CONECT records
        if not line.startswith("CONECT"):
            f_out.write(line)

print(f"Cleaned receptor saved as {output_file}")