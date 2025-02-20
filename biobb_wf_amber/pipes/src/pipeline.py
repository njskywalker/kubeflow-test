import os
import kfp.dsl as dsl
# from kfp import Client
from kfp.compiler import Compiler

from stages.fetch import fetch_pdb_protein
from stages.prepare import prep_pdb_for_amber, prep_amber_topology, prep_amber_to_pdb
from stages.simulate import simulate

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
        input_path=topology_task.outputs['output_path'],
        output_traj_filename='/sander.h_min.x',
        output_rst_filename='/sander.h_min.rst',
        output_log_filename='/sander.h_min.log',
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
        input_path=topology_task.outputs['output_path'],
        output_traj_filename='/sander.n_min.x',
        output_rst_filename='/sander.n_min.rst',
        output_log_filename='/sander.n_min.log',
    )

    amber_to_pdb_task = prep_amber_to_pdb(
        input_topology_path=topology_task.outputs['output_path'],
        input_minimization_path=system_minimization_task.outputs['output_path'],
    )
    # parameterize_task = parameterize_molecule_op(prepare_task.output)
    # gen_box_task = generate_simulation_box_op(parameterize_task.output)
    # minimize_task = energy_minimization_op(gen_box_task.output)


if __name__ == '__main__':
    """Only compiles the pipeline. To run, look at `execute.py`."""

    Compiler().compile(pipeline, '/src/outputs/pipeline.yaml')
