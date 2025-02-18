import os
import kfp.dsl as dsl
# from kfp import Client
from kfp.compiler import Compiler

from stages.fetch import fetch_pdb_protein

@dsl.pipeline(
    name='Isomorphic MD Pipeline',
    description='General pipeline for MD setup given PDB ID.'
)
def pipeline(pdb_code: str):
    fetch_pdb_protein_task = fetch_pdb_protein(pdb_code=pdb_code)
    # fetch_task = fetch_molecular_data_op(molecule_id)
    # prepare_task = prepare_structure_op(fetch_task.output)
    # parameterize_task = parameterize_molecule_op(prepare_task.output)
    # gen_box_task = generate_simulation_box_op(parameterize_task.output)
    # minimize_task = energy_minimization_op(gen_box_task.output)


if __name__ == '__main__':
    """Only compiles the pipeline. To run, look at `execute.py`."""

    Compiler().compile(pipeline, 'pipeline.yaml')
