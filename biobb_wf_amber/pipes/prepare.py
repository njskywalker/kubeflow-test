"""PDB -> AMBER preparation pipe stage."""

from typing import Dict

def prep_pdb_for_amber(downloaded_pdb: str) -> str:
    """Prepares (cleans) fetched PDB files for AMBER analysis.

    For instance, relabels PDB codes for cysteines forming
    disulphide bridges for proper AMBER recognition.
    
    Returns path to processed PDB file."""

    # Import module
    from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run

    # Create prop dict and inputs/outputs
    output_pdb4amber_path = 'structure.pdb4amber.pdb'

    # Create and launch bb
    pdb4amber_run(input_pdb_path=downloaded_pdb,
                output_pdb_path=output_pdb4amber_path)
    
    return output_pdb4amber_path

def prep_amber_topology(input_pdb_path: str) -> Dict[str, str]:
    """Builds AMBER topology for a given PDB protein file path.
    
    Returns a dictionary of three filepaths for output PDB, 
    AMBER topology (top) and coordinate (crd) files."""

    # Import module
    from biobb_amber.leap.leap_gen_top import leap_gen_top

    # Create prop dict and inputs/outputs
    output_pdb_path = 'structure.leap.pdb'
    output_top_path = 'structure.leap.top'
    output_crd_path = 'structure.leap.crd'

    prop = {
        "forcefield" : ["protein.ff14SB"]
    }

    # Create and launch bb
    leap_gen_top(input_pdb_path=input_pdb_path,
            output_pdb_path=output_pdb_path,
            output_top_path=output_top_path,
            output_crd_path=output_crd_path,
            properties=prop)
    
    # TODO: Improve legibility - NamedTuple/Dict or basic data class?
    return {
        "output_pb_path": output_pdb_path, 
        "output_top_path": output_top_path, 
        "output_crd_path": output_crd_path
    }

def prep_amber_to_pdb(
    input_top_path: str,
    input_crd_path: str,
) -> str:
    """Converts AMBER protein representation to PDB.
    
    Returns filepath to """

    # Import module
    from biobb_amber.ambpdb.amber_to_pdb import amber_to_pdb

    # Create prop dict and inputs/outputs
    output_ambpdb_path = 'structure.ambpdb.pdb'

    # Create and launch bb
    amber_to_pdb(input_top_path=input_top_path,
                input_crd_path=input_crd_path,
                output_pdb_path=output_ambpdb_path
                )
    
    return output_ambpdb_path