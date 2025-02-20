import os
import kfp.dsl as dsl
# from kfp import Client
from kfp.compiler import Compiler

from stages.fetch import fetch_pdb_protein
from stages.prepare import prep_pdb_for_amber, prep_amber_topology, prep_amber_to_pdb
from stages.simulate import simulate
from stages.solvate import create_water_box, add_ions

@dsl.pipeline(
    name='Test Kubeflow Pipeline',
    description='General pipeline for MD setup given PDB ID.'
)
def pipeline(pdb_code: str):
    """General pipeline for MD setup given PDB ID."""

    # Fetching
    fetch_pdb_protein_task = fetch_pdb_protein(pdb_code=pdb_code)

    # Data preparation
    prepare_task = prep_pdb_for_amber(input_pdb_dir_path=fetch_pdb_protein_task.outputs['output_path'])
    topology_task = prep_amber_topology(
        input_path=prepare_task.outputs['output_amber_dir_path'],
        properties={"forcefield" : ["protein.ff14SB"]}
    )

    # Minimization / simulation steps
    # Brings system closer to "feasible" "realistic" starting point
    hydrogen_minimization_task = simulate(
        properties = {
            'simulation_type' : "min_vacuo",
            "mdin" : { 
                'maxcyc' : 500,
                'ntpr' : 5,
                'ntr' : 1,
                'restraintmask' : '\":*&!@H=\"',
                'restraint_wt' : 50.0
            }
        },
        input_top_path=topology_task.outputs['output_path'],
        input_crd_path=topology_task.outputs['output_path'],
        input_top_filename="structure.leap.top",
        input_crd_filename="structure.leap.crd",
        output_traj_filename='sander.h_min.x',
        output_rst_filename='sander.h_min.rst',
        output_log_filename='sander.h_min.log',
    )
    # TODO: parallelisable analysis steps (to create graphs, pngs etc.)

    # TODO: This is exactly the same task as above, is this just
    # a tutorial thing or is this a mistake? No diffs in props
    system_minimization_task = simulate(
        properties = {
            'simulation_type' : "min_vacuo",
            "mdin" : { 
                'maxcyc' : 500,
                'ntpr' : 5,
                'ntr' : 1,
                'restraintmask' : '\":*&!@H=\"',
                'restraint_wt' : 50.0
            }
            },  
        input_top_path=topology_task.outputs['output_path'],
        input_crd_path=topology_task.outputs['output_path'],
        output_traj_filename='sander.n_min.x',
        output_rst_filename='sander.n_min.rst',
        output_log_filename='sander.n_min.log',
        input_top_filename="structure.leap.top",
        input_crd_filename="structure.leap.crd",
    )

    amber_to_pdb_task = prep_amber_to_pdb(
        input_topology_path=topology_task.outputs['output_path'],
        input_minimization_path=system_minimization_task.outputs['output_path'],
    )
    
    # solvate
    water_box_task = create_water_box(
        properties = {
            "forcefield" : ["protein.ff14SB"],
            "water_type": "TIP3PBOX",
            "distance_to_molecule": "9.0",   
            "box_type": "truncated_octahedron"
        },
        input_path=amber_to_pdb_task.outputs['output_path'],
        output_solv_pdb_filename='structure.solv.pdb',
        output_solv_top_filename='structure.solv.parmtop',
        output_solv_crd_filename='structure.solv.crd',
    )

    # ionate
    add_ions_task = add_ions(
        properties={
            "forcefield" : ["protein.ff14SB"],
            "neutralise" : True,
            "positive_ions_type": "Na+",
            "negative_ions_type": "Cl-",
            "ionic_concentration" : 150, # 150mM
            "box_type": "truncated_octahedron"
        },
        input_path=water_box_task.outputs['output_path'],
        output_ions_pdb_filename='structure.ions.pdb',
        output_ions_top_filename='structure.ions.parmtop',
        output_ions_crd_filename='structure.ions.crd',
    )

    # minimize
    steepest_descent_task = simulate(
        properties = {
            "simulation_type" : "minimization",
            "mdin" : { 
                'maxcyc' : 300, # Reducing the number of minimization steps for the sake of time
                'ntr' : 1,      # Overwritting restrain parameter
                'restraintmask' : '\"!:WAT,Cl-,Na+\"',      # Restraining solute
                'restraint_wt' : 50.0                       # With a force constant of 50 Kcal/mol*A2
            }
        },
        input_top_path=add_ions_task.outputs['output_path'],
        input_crd_path=add_ions_task.outputs['output_path'],
        output_traj_filename='sander.min.x',
        output_rst_filename='sander.min.rst',
        output_log_filename='sander.min.log',
        input_top_filename="structure.ions.parmtop",
        input_crd_filename="structure.ions.crd",
    )

    # heat
    heating_task = simulate(
        properties={
            "simulation_type" : "heat",
            "mdin" : { 
                'nstlim' : 2500, # Reducing the number of steps for the sake of time (5ps)
                'ntr' : 1,       # Overwritting restrain parameter
                'restraintmask' : '\"!:WAT,Cl-,Na+\"',      # Restraining solute
                'restraint_wt' : 10.0                       # With a force constant of 10 Kcal/mol*A2
            }
        },
        output_traj_filename='sander.heat.netcdf',
        output_rst_filename='sander.heat.rst',
        output_log_filename='sander.heat.log',
        input_top_path=add_ions_task.outputs['output_path'],
        input_top_filename="structure.ions.parmtop",
        input_crd_path=steepest_descent_task.outputs['output_path'],
        input_crd_filename="sander.min.rst",
    )

    # nvt
    nvt_simulation_task = simulate(
        properties={
            "simulation_type" : 'nvt',
            "mdin" : { 
                'nstlim' : 500, # Reducing the number of steps for the sake of time (1ps)
                'ntr' : 1,      # Overwritting restrain parameter
                'restraintmask' : '\"!:WAT,Cl-,Na+ & !@H=\"',      # Restraining solute heavy atoms
                'restraint_wt' : 5.0                               # With a force constant of 5 Kcal/mol*A2
            }
        },
        output_traj_filename='sander.nvt.netcdf',
        output_rst_filename='sander.nvt.rst',
        output_log_filename='sander.nvt.log',
        input_top_path=add_ions_task.outputs['output_path'],
        input_top_filename="structure.ions.parmtop",
        input_crd_path=heating_task.outputs['output_path'],
        input_crd_filename="sander.heat.rst",
    )

    # npt
    npt_simulation_task = simulate(
        properties={
            "simulation_type" : 'npt',
            "mdin" : { 
                'nstlim' : 500, # Reducing the number of steps for the sake of time (1ps)
                'ntr' : 1,      # Overwritting restrain parameter
                'restraintmask' : '\"!:WAT,Cl-,Na+ & !@H=\"',      # Restraining solute heavy atoms
                'restraint_wt' : 2.5                               # With a force constant of 2.5 Kcal/mol*A2
            }
        },
        output_traj_filename='sander.npt.netcdf',
        output_rst_filename='sander.npt.rst',
        output_log_filename='sander.npt.log',
        input_top_path=add_ions_task.outputs['output_path'],
        input_top_filename="structure.ions.parmtop",
        input_crd_path=nvt_simulation_task.outputs['output_path'],
        input_crd_filename="sander.nvt.rst",
    )

    # simulate (Free MD)
    free_simulation_task = simulate(
        properties={
            "simulation_type" : 'free',
            "mdin" : { 
                'nstlim' : 2500, # Reducing the number of steps for the sake of time (5ps)
                'ntwx' : 500  # Print coords to trajectory every 500 steps (1 ps)
            }
        },
        output_traj_filename='sander.free.netcdf',
        output_rst_filename='sander.free.rst',
        output_log_filename='sander.free.log',
        input_top_path=add_ions_task.outputs['output_path'],
        input_top_filename="structure.ions.parmtop",
        input_crd_path=npt_simulation_task.outputs['output_path'],
        input_crd_filename="sander.npt.rst",
    )



if __name__ == '__main__':
    """Only compiles the pipeline. To run, look at `execute.py`."""

    Compiler().compile(pipeline, '/src/outputs/pipeline.yaml')
