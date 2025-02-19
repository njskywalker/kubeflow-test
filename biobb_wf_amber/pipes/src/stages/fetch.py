from kfp import dsl
from src.common.artifacts import PDBFile
import subprocess

@dsl.component(
        packages_to_install=['biobb_io'],
        base_image='pipeline:latest',
)
def fetch_pdb_protein(pdb_code: str, output_path: dsl.Output[PDBFile]) -> None:
    """Fetches a PDB protein using its code.
    
    Returns path to downloaded PDB file."""
    
    import os
    from biobb_io.api.pdb import pdb

    prop = {
        'pdb_code': pdb_code
    }

    # Create and launch bb
    pdb(output_pdb_path=output_path,
        properties=prop)
