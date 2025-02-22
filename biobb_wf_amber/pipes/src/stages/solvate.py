"""Solvation pipeline step logic."""

from kfp import dsl

@dsl.component(
        # packages_to_install=['biobb_amber'],
        base_image="amber_bio:latest"
)
def create_water_box(
        output_solv_pdb_filename: str,
        output_solv_top_filename: str,
        output_solv_crd_filename: str,
        properties: dict,
        input_path: dsl.InputPath('Directory'),
        output_path: dsl.OutputPath('Directory'),
):
    """Creates a unit cell around protein (PBD input filepath)
    and solvates with water."""

    # Import modules
    import os
    from biobb_amber.leap.leap_solvate import leap_solvate

    if os.getenv("AMBERHOME") is None:
        os.environ["AMBERHOME"] = "/opt/conda"

    # Inputs and outputs
    input_pdb_path = input_path + "/" + 'structure.ambpdb.pdb'
    os.makedirs(output_path, exist_ok=True)
    output_solv_pdb_path = output_path + "/" + output_solv_pdb_filename
    output_solv_top_path = output_path + "/" + output_solv_top_filename
    output_solv_crd_path = output_path + "/" + output_solv_crd_filename

    {
        "forcefield" : ["protein.ff14SB"],
        "water_type": "TIP3PBOX",
        "distance_to_molecule": "9.0",   
        "box_type": "truncated_octahedron"
    }

    # Create and launch bb
    leap_solvate(
        input_pdb_path=input_pdb_path,
        output_pdb_path=output_solv_pdb_path,
        output_top_path=output_solv_top_path,
        output_crd_path=output_solv_crd_path,
        properties=properties
    )


@dsl.component(
        packages_to_install=['biobb_amber'],
        base_image="amber_bio:latest"
)
def add_ions(
        properties: dict,
        output_ions_pdb_filename: str,
        output_ions_top_filename: str,
        output_ions_crd_filename: str,
        input_path: dsl.InputPath('Directory'),
        output_path: dsl.OutputPath('Directory'),
):
    """Adds ions in solvated water box."""

    # Import modules
    import os
    from biobb_amber.leap.leap_add_ions import leap_add_ions

    if os.getenv("AMBERHOME") is None:
        os.environ["AMBERHOME"] = "/opt/conda"

    # Create prop dict and inputs/outputs
    input_pdb_path = input_path + "/" + "structure.solv.pdb"
    os.makedirs(output_path, exist_ok=True)
    output_ions_pdb_path = output_path + "/" + output_ions_pdb_filename
    output_ions_top_path = output_path + "/" + output_ions_top_filename
    output_ions_crd_path = output_path + "/" + output_ions_crd_filename

    # Create and launch bb
    leap_add_ions(input_pdb_path=input_pdb_path,
            output_pdb_path=output_ions_pdb_path,
            output_top_path=output_ions_top_path,
            output_crd_path=output_ions_crd_path,
            properties=properties)
