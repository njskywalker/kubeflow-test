from kfp import dsl

@dsl.component(
        packages_to_install=['biobb_io']
)
def fetch_pdb_protein(pdb_code: str) -> str:
    """Fetches a PDB protein using its code.
    
    Returns path to downloaded PDB file."""
    
    # Import module
    from biobb_io.api.pdb import pdb

    # Create properties dict and inputs/outputs
    downloaded_pdb_filename = '/src/outputs/' + pdb_code + '.pdb'

    prop = {
        'pdb_code': pdb_code
    }

    # Create and launch bb
    pdb(output_pdb_path=downloaded_pdb_filename,
        properties=prop)

    return downloaded_pdb_filename
