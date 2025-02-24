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
        # packages_to_install=['biobb_amber'],
        # base_image="quay.io/biocontainers/biobb_amber:5.0.4--pyhdfd78af_0",
        # above also doesn't work, but combining with own image
        # installs `tleap` ctl which is needed for `leap_gen_top`
        base_image="amber_bio:latest"
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
    import subprocess

    # Workaround because image is buggy
    # Not sure how it loads conda venv (or doesn't)
    # Perhaps custom entrypoint used? Flushing env?
    # Need more debugging
    if os.getenv("AMBERHOME") is None:
        os.environ["AMBERHOME"] = "/opt/conda"

    # Paths
    os.makedirs(output_path, exist_ok=True)
    input_pdb_path = input_path + '/structure.pdb4amber.pdb'
    output_pdb_path = output_path + '/structure.leap.pdb'
    output_top_path = output_path + '/structure.leap.top'
    output_crd_path = output_path + '/structure.leap.crd'

    # Create and launch bb
    leap_gen_top(
        input_pdb_path=input_pdb_path, 
        output_pdb_path=output_pdb_path, 
        output_top_path=output_top_path, 
        output_crd_path=output_crd_path, 
        properties=properties
    )
    

@dsl.component(
        # packages_to_install=['biobb_amber'],
        base_image="amber_bio:latest"
)
def prep_amber_to_pdb(
    input_topology_path: dsl.InputPath('Directory'),
    input_minimization_path: dsl.InputPath('Directory'),
    output_path: dsl.OutputPath('Directory'),
):
    """Converts AMBER protein representation to PDB.
    
    Returns filepath to """

    # Import modules
    import os
    from biobb_amber.ambpdb.amber_to_pdb import amber_to_pdb

    # Create input/output paths
    input_top_path = input_topology_path + "/structure.leap.top"
    input_crd_path = input_minimization_path + "/sander.n_min.rst"
    os.makedirs(output_path, exist_ok=True)
    output_ambpdb_path = output_path + '/structure.ambpdb.pdb'

    # Create and launch bb
    amber_to_pdb(
        input_top_path=input_top_path,
        input_crd_path=input_crd_path,
        output_pdb_path=output_ambpdb_path
    )
