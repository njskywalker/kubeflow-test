"""PDB -> AMBER preparation pipe stage."""

from typing import Dict, Any
from kfp import dsl


@dsl.component(
        # packages_to_install=['biobb_amber'],  # works even tho pip install
        # but base image faster
        base_image="quay.io/biocontainers/biobb_amber:5.0.4--pyhdfd78af_0"
)
def prep_pdb_for_amber(
    input_pdb_dir_path: dsl.InputPath('Directory'),
    output_amber_dir_path: dsl.OutputPath('Directory')
) -> None:
    """Prepares (cleans) fetched PDB files for AMBER analysis.

    For instance, relabels PDB codes for cysteines forming
    disulphide bridges for proper AMBER recognition."""

    # Import modules
    import os
    from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run

    # Create input/output paths
    input_pdb_path = input_pdb_dir_path + '/protein.pdb'
    os.makedirs(output_amber_dir_path, exist_ok=True)
    output_pdb_path = output_amber_dir_path + '/structure.pdb4amber.pdb'

    # Create and launch bb
    pdb4amber_run(input_pdb_path=input_pdb_path,
                output_pdb_path=output_pdb_path)



@dsl.component(
        packages_to_install=['biobb_amber'],
        # base_image="quay.io/biocontainers/biobb_amber:5.0.4--pyhdfd78af_0",
        # above also doesn't work
        base_image="nebjovanovic/amber_bio:latest"
)
def prep_amber_topology(
    properties: Dict[str, Any],
    input_path: dsl.InputPath('Directory'),
    output_path: dsl.OutputPath('Directory'),
) -> None:
    """Builds AMBER topology for a given PDB protein file path.
    
    Returns a dictionary of three filepaths for output PDB, 
    AMBER topology (top) and coordinate (crd) files."""

    # Import module
    import os
    from biobb_amber.leap.leap_gen_top import leap_gen_top

    print("Printing inside container func")

    # Paths
    os.makedirs(output_path, exist_ok=True)
    input_pdb_path = '/structure.pdb4amber.pdb'
    output_pdb_path = '/structure.leap.pdb'
    output_top_path = '/structure.leap.top'
    output_crd_path = '/structure.leap.crd'
    properties={"forcefield" : ["protein.ff14SB"]}

    # Create and launch bb
    leap_gen_top(input_pdb_path="/proteinamb.pdb", output_pdb_path=output_pdb_path, output_top_path=output_top_path, output_crd_path=output_crd_path, properties={"forcefield" : ["protein.ff14SB"]})
    
    # # TODO: Improve legibility - NamedTuple/Dict or basic data class?
    # return {
    #     "output_pb_path": output_pdb_path, 
    #     "output_top_path": output_top_path, 
    #     "output_crd_path": output_crd_path
    # }

# def prep_amber_to_pdb(
#     input_top_path: str,
#     input_crd_path: str,
# ) -> str:
#     """Converts AMBER protein representation to PDB.
    
#     Returns filepath to """

#     # Import module
#     from biobb_amber.ambpdb.amber_to_pdb import amber_to_pdb

#     # Create prop dict and inputs/outputs
#     output_ambpdb_path = 'structure.ambpdb.pdb'

#     # Create and launch bb
#     amber_to_pdb(input_top_path=input_top_path,
#                 input_crd_path=input_crd_path,
#                 output_pdb_path=output_ambpdb_path
#                 )
    
#     return output_ambpdb_path
