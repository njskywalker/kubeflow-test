from kfp import dsl
from pipes.artifacts import PDBFile

@dsl.component(
        packages_to_install=['biobb_io']
)
def fetch_pdb_protein(pdb_code: str, output_path: dsl.Output[PDBFile]) -> None:
    """Fetches a PDB protein using its code.
    
    Returns path to downloaded PDB file."""
    
    import os
    from biobb_io.api.pdb import pdb

    prop = {
        'pdb_code': pdb_code
    }

    # if not output_path.endswith('.pdb'):
    #     output_path = output_path + '.pdb'

    # Create and launch bb
    pdb(output_pdb_path=output_path,
        properties=prop)
