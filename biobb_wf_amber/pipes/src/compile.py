import os
import kfp.dsl as dsl
# from kfp import Client
from kfp.compiler import Compiler


from stages.prepare import prep_pdb_for_amber

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


@dsl.pipeline(
    name='Test Kubeflow Pipeline',
    description='General pipeline for MD setup given PDB ID.'
)
def pipeline(pdb_code: str):
    fetch_pdb_protein_task = fetch_pdb_protein(pdb_code=pdb_code)
    prepare_task = prep_pdb_for_amber(input_pdb_dir_path=fetch_pdb_protein_task.outputs['output_path'])
    # parameterize_task = parameterize_molecule_op(prepare_task.output)
    # gen_box_task = generate_simulation_box_op(parameterize_task.output)
    # minimize_task = energy_minimization_op(gen_box_task.output)


if __name__ == '__main__':
    """Only compiles the pipeline. To run, look at `execute.py`."""

    Compiler().compile(pipeline, '/src/outputs/pipeline.yaml')
