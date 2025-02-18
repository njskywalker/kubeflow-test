"""Energy minimisation pipe functions."""

def minimize_hydrogen(
    input_top_path: str,
    input_crd_path: str
) -> Dict[str, str]:
    """Optimize hydrogen positions of given AMBER protein representation.
    
    Takes in AMBER topology and coordinate (top and crd) filepaths.
    Returns dictionary with output filepaths.
    """

    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_h_min_traj_path = 'sander.h_min.x'
    output_h_min_rst_path = 'sander.h_min.rst'
    output_h_min_log_path = 'sander.h_min.log'

    prop = {
        'simulation_type' : "min_vacuo",
        "mdin" : { 
            'maxcyc' : 500,
            'ntpr' : 5,
            'ntr' : 1,
            'restraintmask' : '\":*&!@H=\"',
            'restraint_wt' : 50.0
        }
    }

    # Create and launch bb
    sander_mdrun(input_top_path=input_top_path,
                input_crd_path=input_crd_path,
                input_ref_path=input_crd_path,
                output_traj_path=output_h_min_traj_path,
                output_rst_path=output_h_min_rst_path,
                output_log_path=output_h_min_log_path,
                properties=prop)

    return {
        "output_traj_path": output_h_min_traj_path,
        "output_rst_path": output_h_min_rst_path,
        "output_log_path": output_h_min_log_path,
    }