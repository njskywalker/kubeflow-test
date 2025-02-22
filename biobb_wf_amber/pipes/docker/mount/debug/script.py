# Import module
import os
from biobb_amber.leap.leap_gen_top import leap_gen_top

print("In the latest")

# Paths
input_pdb_path = 'structure.pdb4amber.pdb'
output_pdb_path = 'structure.leap.pdb'
output_top_path = 'structure.leap.top'
output_crd_path = 'structure.leap.crd'

print(os.getenv("AMBERHOME"))

# Create and launch bb
leap_gen_top(
    input_pdb_path=input_pdb_path, 
    output_pdb_path=output_pdb_path, 
    output_top_path=output_top_path, 
    output_crd_path=output_crd_path, 
    properties={"forcefield" : ["protein.ff14SB"]}
)
