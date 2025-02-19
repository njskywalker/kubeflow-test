"""Solvation pipeline step logic."""

def create_water_box(
        input_pdb_path: str,
) -> str:
    """Creates a unit cell around protein (PBD input filepath)
    and solvates with water."""

    # Import module
    from biobb_amber.leap.leap_solvate import leap_solvate

    # Create prop dict and inputs/outputs
    output_solv_pdb_path = 'structure.solv.pdb'
    output_solv_top_path = 'structure.solv.parmtop'
    output_solv_crd_path = 'structure.solv.crd'

    prop = {
        "forcefield" : ["protein.ff14SB"],
        "water_type": "TIP3PBOX",
        "distance_to_molecule": "9.0",   
        "box_type": "truncated_octahedron"
    }

    # Create and launch bb
    leap_solvate(input_pdb_path=input_pdb_path,
            output_pdb_path=output_solv_pdb_path,
            output_top_path=output_solv_top_path,
            output_crd_path=output_solv_crd_path,
            properties=prop)
    
    return {
        "output_solv_pdb_path": output_solv_pdb_path,
        "output_solv_top_path": output_solv_top_path,
        "output_solv_crd_path": output_solv_crd_path,
    }

def add_ions(
        pdb_path: str,
        positive_ions_type: str = "Na+",
        negative_ions_type: str = "Cl-",
) -> str:
    """Adds ions in solvated water box."""

    # Import module
    from biobb_amber.leap.leap_add_ions import leap_add_ions

    # Create prop dict and inputs/outputs
    output_ions_pdb_path = 'structure.ions.pdb'
    output_ions_top_path = 'structure.ions.parmtop'
    output_ions_crd_path = 'structure.ions.crd'

    prop = {
        "forcefield" : ["protein.ff14SB"],
        "neutralise" : True,
        "positive_ions_type": positive_ions_type, # "Na+",
        "negative_ions_type": negative_ions_type, # "Cl-",
        "ionic_concentration" : 150, # 150mM
        "box_type": "truncated_octahedron"
    }

    # Create and launch bb
    leap_add_ions(input_pdb_path=pdb_path,
            output_pdb_path=output_ions_pdb_path,
            output_top_path=output_ions_top_path,
            output_crd_path=output_ions_crd_path,
            properties=prop)
    
    return {
        "output_ions_pdb_path": output_ions_pdb_path, 
        "output_ions_top_path": output_ions_top_path, 
        "output_ions_crd_path": output_ions_crd_path, 
    }