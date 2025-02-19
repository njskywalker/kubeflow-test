from kfp import dsl

@dsl.component(
        packages_to_install=['biobb_io'],
)
def fetch_pdb_protein(pdb_code: str, output_path: dsl.OutputPath('Directory')) -> None:
    """Fetches a PDB protein using its code.
    
    Returns path to downloaded PDB file."""

    import os
    from biobb_io.api.pdb import pdb

    prop = {
        'pdb_code': pdb_code
    }

    os.makedirs(output_path, exist_ok=True)
    output_pdb_path = output_path + '/protein.pdb'

    # Create and launch bb
    pdb(output_pdb_path=output_pdb_path,
        properties=prop)
