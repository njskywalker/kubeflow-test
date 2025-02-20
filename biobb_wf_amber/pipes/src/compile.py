import os
import kfp.dsl as dsl
# from kfp import Client
from kfp.compiler import Compiler

from stages.fetch import fetch_pdb_protein
from stages.prepare import prep_pdb_for_amber, prep_amber_topology

@dsl.pipeline(
    name='Test Kubeflow Pipeline',
    description='General pipeline for MD setup given PDB ID.'
)
def pipeline(pdb_code: str):
    """General pipeline for MD setup given PDB ID."""

    # Fetching
    fetch_pdb_protein_task = fetch_pdb_protein(pdb_code=pdb_code)

    # Preparation
    prepare_task = prep_pdb_for_amber(input_pdb_dir_path=fetch_pdb_protein_task.outputs['output_path'])
    topology_task = prep_amber_topology(
        input_path=prepare_task.outputs['output_amber_dir_path'],
        properties={"forcefield" : ["protein.ff14SB"]}
    )

    # parameterize_task = parameterize_molecule_op(prepare_task.output)
    # gen_box_task = generate_simulation_box_op(parameterize_task.output)
    # minimize_task = energy_minimization_op(gen_box_task.output)


if __name__ == '__main__':
    """Only compiles the pipeline. To run, look at `execute.py`."""

    Compiler().compile(pipeline, '/src/outputs/pipeline.yaml')
