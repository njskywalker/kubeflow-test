"""Energy minimisation pipe functions."""

from typing import Dict, Any
from kfp import dsl

# def _simulate(
#     input_top_path_complete: str,
#     input_crd_path_complete: str,
#     output_path: str,
#     output_traj_filename: str,
#     output_rst_filename: str,
#     output_log_filename: str,
#     properties: Dict[str, Any]
# ):
    


@dsl.component(
        # packages_to_install=['biobb_amber'],
        base_image="amber_bio:latest"
)
def simulate_one_input(
    properties: Dict[str, Any],
    output_traj_filename: str,
    output_rst_filename: str,
    output_log_filename: str,
    input_top_path: dsl.InputPath('Directory'),
    input_top_filename: str,
    # input_crd_path: dsl.InputPath('Directory'),
    input_crd_filename: str,
    output_path: dsl.OutputPath('Directory'),
):
    """
    Run a simulation with given parameters with AmberTools Sander.

    Accepts only one input path, as a Kubeflow Go component throws an error 
    if the same InputPath is passed to two different function InputPath parameters.
    Example:

    # KFP driver: driver.Container(
    #     pipelineName=test-kubeflow-pipeline, runID=24e1e2c6-e639-4f87-97ab-953608807c2e, 
    #     task="simulate", component="comp-simulate", dagExecutionID=44, componentSpec) 
    #     failed: rpc error: code = AlreadyExists desc = Given event already exists: 
    #     artifact_id: 28
    # execution_id: 48
    # path {
    #     steps {
    #         key: "input_top_path_x"
    #     }
    # }

    Inputs:
    - properties: Dict containing simulation parameters.
    - output_traj_filename: Desired name of trajectory file. (NB May not be created)
    - output_rst_filename: Desired name of restart file.
    - output_log_filename: Desired name of log file.
    - input_path: Path to input files from previous stage.

    Creates folder with outputs.
    """

    # Modify inputs/outputs
    input_top_path_complete = input_top_path + "/" + input_top_filename # "structure.leap.top"
    input_crd_path_complete = input_top_path + "/" + input_crd_filename # "structure.leap.crd"

    # Import shenanigans
    # from simulate import _simulate

    # Import module
    import os
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    os.makedirs(output_path, exist_ok=True)

    # TODO: look into using os.path.join - but struggles with slashes?
    # or just add a `/` here...
    output_h_min_traj_path = output_path + "/" + output_traj_filename
    output_h_min_rst_path = output_path + "/" + output_rst_filename
    output_h_min_log_path = output_path + "/" + output_log_filename


    # Create and launch bb
    sander_mdrun(
        input_top_path=input_top_path_complete,
        input_crd_path=input_crd_path_complete,
        input_ref_path=input_crd_path_complete,
        output_traj_path=output_h_min_traj_path,
        output_rst_path=output_h_min_rst_path,
        output_log_path=output_h_min_log_path,
        properties=properties
    )


@dsl.component(
        # packages_to_install=['biobb_amber'],
        base_image="amber_bio:latest"
)
def simulate_two_inputs(
    properties: Dict[str, Any],
    output_traj_filename: str,
    output_rst_filename: str,
    output_log_filename: str,
    input_top_path: dsl.InputPath('Directory'),
    input_top_filename: str,
    input_crd_path: dsl.InputPath('Directory'),
    input_crd_filename: str,
    output_path: dsl.OutputPath('Directory'),
):
    """
    The same as simulate_one_input, but with two input paths which
    must be different.

    Required as Kubeflow throws an error if the same InputPath is passed
    to two different function InputPath parameters.
    """

    # Modify inputs/outputs
    input_top_path_complete = input_top_path + "/" + input_top_filename # "structure.leap.top"
    input_crd_path_complete = input_crd_path + "/" + input_crd_filename # "structure.leap.crd"

    # Import shenanigans
    # from simulate import _simulate

    # Import module
    import os
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    os.makedirs(output_path, exist_ok=True)

    # TODO: look into using os.path.join - but struggles with slashes?
    # or just add a `/` here...
    output_h_min_traj_path = output_path + "/" + output_traj_filename
    output_h_min_rst_path = output_path + "/" + output_rst_filename
    output_h_min_log_path = output_path + "/" + output_log_filename


    # Create and launch bb
    sander_mdrun(
        input_top_path=input_top_path_complete,
        input_crd_path=input_crd_path_complete,
        input_ref_path=input_crd_path_complete,
        output_traj_path=output_h_min_traj_path,
        output_rst_path=output_h_min_rst_path,
        output_log_path=output_h_min_log_path,
        properties=properties
    )
