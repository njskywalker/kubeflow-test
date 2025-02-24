import kfp.dsl as dsl
from kfp.compiler import Compiler
from kfp import kubernetes

from stages.fetch import fetch_pdb_protein
from stages.prepare import prep_pdb_for_amber, prep_amber_topology, prep_amber_to_pdb
from stages.simulate import simulate_one_input, simulate_two_inputs
from stages.solvate import create_water_box, add_ions
from analysis.scatter import plot_simulation_scatter


@dsl.pipeline(
    name='PDB -> AMBER protein processing pipeline.',
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
        properties={"forcefield" : ["protein.ff14SB"]},
    )
    kubernetes.set_image_pull_policy(task=topology_task, policy="Never")

    # Minimization / simulation steps
    # Brings system closer to "feasible" "realistic" starting point
    hydrogen_minimization_task = simulate_one_input(
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
        # input_crd_path=topology_task.outputs['output_path'],
        input_top_filename="structure.leap.top",
        input_crd_filename="structure.leap.crd",
        output_traj_filename='sander.h_min.x',
        output_rst_filename='sander.h_min.rst',
        output_log_filename='sander.h_min.log',
    )
    kubernetes.set_image_pull_policy(task=hydrogen_minimization_task, policy="Never")

    process_hydrogen_minimization_task = plot_simulation_scatter(
        properties={
            "terms" : ['ENERGY']
        },
        input_path=hydrogen_minimization_task.outputs['output_path'],
        input_log_filename='sander.h_min.log',
        output_dat_filename='sander.h_min.energy.dat',
        plot_title='Energy Minimization',
        x_axis_title='Energy Minimization Step',
        y_axis_title='Potential Energy (Kcal/mol)',
        plot_height=600,
        process_type='minout',
    )
    kubernetes.set_image_pull_policy(task=process_hydrogen_minimization_task, policy="Never")

    # TODO: This is exactly the same task as above, is this just
    # a tutorial thing or is this a mistake? No diffs in props
    system_minimization_task = simulate_one_input(
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
        # input_crd_path=topology_task.outputs['output_path'],
        output_traj_filename='sander.n_min.x',
        output_rst_filename='sander.n_min.rst',
        output_log_filename='sander.n_min.log',
        input_top_filename="structure.leap.top",
        input_crd_filename="structure.leap.crd",
    )
    kubernetes.set_image_pull_policy(task=system_minimization_task, policy="Never")
    process_system_minimization_task = plot_simulation_scatter(
        properties={
            "terms" : ['ENERGY']
        },
        input_path=system_minimization_task.outputs['output_path'],
        input_log_filename='sander.n_min.log',
        output_dat_filename='sander.n_min.energy.dat',

        # would change the following to be more specific
        # but not touching how notebook owner wants it!
        plot_title='Energy Minimization',
        x_axis_title='Energy Minimization Step',
        y_axis_title='Potential Energy (Kcal/mol)',
        plot_height=600,
        process_type='minout',
    )
    kubernetes.set_image_pull_policy(task=process_system_minimization_task, policy="Never")

    amber_to_pdb_task = prep_amber_to_pdb(
        input_topology_path=topology_task.outputs['output_path'],
        input_minimization_path=system_minimization_task.outputs['output_path'],
    )
    kubernetes.set_image_pull_policy(task=amber_to_pdb_task, policy="Never")
    
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
    kubernetes.set_image_pull_policy(task=water_box_task, policy="Never")

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
    kubernetes.set_image_pull_policy(task=add_ions_task, policy="Never")

    # minimize
    steepest_descent_task = simulate_one_input(
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
        # input_crd_path=add_ions_task.outputs['output_path'],
        output_traj_filename='sander.min.x',
        output_rst_filename='sander.min.rst',
        output_log_filename='sander.min.log',
        input_top_filename="structure.ions.parmtop",
        input_crd_filename="structure.ions.crd",
    )
    kubernetes.set_image_pull_policy(task=steepest_descent_task, policy="Never")
    process_steepest_descent_task = plot_simulation_scatter(
        properties={
            "terms" : ['ENERGY']
        },
        input_path=steepest_descent_task.outputs['output_path'],
        input_log_filename='sander.min.log',
        output_dat_filename='sander.min.energy.dat',

        # would change the following to be more specific
        # but not touching how notebook owner wants it!
        plot_title='Energy Minimization',
        x_axis_title='Energy Minimization Step',
        y_axis_title='Potential Energy (Kcal/mol)',
        plot_height=600,
        process_type='minout',
    )
    kubernetes.set_image_pull_policy(task=process_steepest_descent_task, policy="Never")

    # heat
    heating_task = simulate_two_inputs(
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
    kubernetes.set_image_pull_policy(task=heating_task, policy="Never")
    process_heating_task = plot_simulation_scatter(
        properties={
            "terms" : ['TEMP']
        },
        input_path=heating_task.outputs['output_path'],
        input_log_filename='sander.heat.log',
        output_dat_filename='sander.md.temp.dat',
        plot_title='Heating process',
        x_axis_title='Heating Step (ps)',
        y_axis_title='Temperature (K)',
        plot_height=600,
        process_type='mdout',
    )
    kubernetes.set_image_pull_policy(task=process_heating_task, policy="Never")

    # nvt
    nvt_simulation_task = simulate_two_inputs(
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
    kubernetes.set_image_pull_policy(task=nvt_simulation_task, policy="Never")
    process_nvt_simulation_task = plot_simulation_scatter(
        properties={
            "terms" : ['TEMP']
        },
        input_path=nvt_simulation_task.outputs['output_path'],
        input_log_filename='sander.nvt.log',
        output_dat_filename='sander.md.nvt.temp.dat',
        plot_title='NVT equilibration',
        x_axis_title='Equilibration Step (ps)',
        y_axis_title='Temperature (K)',
        plot_height=600,
        process_type='mdout',
    )
    kubernetes.set_image_pull_policy(task=process_nvt_simulation_task, policy="Never")

    # npt
    npt_simulation_task = simulate_two_inputs(
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
    kubernetes.set_image_pull_policy(task=npt_simulation_task, policy="Never")

    # simulate (Free MD)
    free_simulation_task = simulate_two_inputs(
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
    kubernetes.set_image_pull_policy(task=free_simulation_task, policy="Never")



if __name__ == '__main__':
    """Only compiles the pipeline. To run, look at `execute.py`."""

    Compiler().compile(pipeline, '/src/outputs/pipeline.yaml')
