"""Energy minimisation pipe functions."""

from typing import Dict, Union

def minimize_hydrogen(
    input_top_path: str,
    input_crd_path: str,
    prop: Dict[str, Union[str, Dict]]
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

    # prop = {
    #     'simulation_type' : "min_vacuo",
    #     "mdin" : { 
    #         'maxcyc' : 500,
    #         'ntpr' : 5,
    #         'ntr' : 1,
    #         'restraintmask' : '\":*&!@H=\"',
    #         'restraint_wt' : 50.0
    #     }
    # }

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


def minimize_system(
        input_top_path: str,
        input_h_min_rst_path: str,
        prop: Dict[str, Union[str, Dict]]
):
    # TODO: Protein-heavy atom minimization seems the same as 
    # hydrogen atom minimization? Speak to data scientist
    
    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_n_min_traj_path = 'sander.n_min.x'
    output_n_min_rst_path = 'sander.n_min.rst'
    output_n_min_log_path = 'sander.n_min.log'

    # prop = {
    #     'simulation_type' : "min_vacuo",
    #     "mdin" : { 
    #         'maxcyc' : 500,
    #         'ntpr' : 5,
    #         'ntr' : 1,
    #         'restraintmask' : '\":*&!@H=\"',
    #         'restraint_wt' : 50.0
    #     }
    # }

    # Create and launch bb
    sander_mdrun(input_top_path=input_top_path,  # output_top_path
                input_crd_path=input_h_min_rst_path,  # output_h_min_rst_path
                input_ref_path=input_h_min_rst_path,
                output_traj_path=output_n_min_traj_path,
                output_rst_path=output_n_min_rst_path,
                output_log_path=output_n_min_log_path,
                properties=prop)


def minimize_steepest_descent(
        input_top_path: str,
        input_crd_path: str,
):
    """Minimize entire system using steepest descent algorithm.
    
    Input should be solvated, ion-filled water box with protein.
    """

    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_min_traj_path = 'sander.min.x'
    output_min_rst_path = 'sander.min.rst'
    output_min_log_path = 'sander.min.log'

    prop = {
        "simulation_type" : "minimization",
        "mdin" : { 
            'maxcyc' : 300, # Reducing the number of minimization steps for the sake of time
            'ntr' : 1,      # Overwritting restrain parameter
            'restraintmask' : '\"!:WAT,Cl-,Na+\"',      # Restraining solute
            'restraint_wt' : 50.0                       # With a force constant of 50 Kcal/mol*A2
        }
    }

    # Create and launch bb
    sander_mdrun(input_top_path=input_top_path,
                input_crd_path=input_crd_path,
                input_ref_path=input_crd_path,
                output_traj_path=output_min_traj_path,
                output_rst_path=output_min_rst_path,
                output_log_path=output_min_log_path,
                properties=prop)
    
    return {
        "output_min_traj_path": output_min_traj_path,
        "output_min_rst_path": output_min_rst_path,
        "output_min_log_path": output_min_log_path,
    }


def heat_system(
        input_top_path: str,
        input_crd_path: str,
):
    """Heats up given system."""

    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_heat_traj_path = 'sander.heat.netcdf'
    output_heat_rst_path = 'sander.heat.rst'
    output_heat_log_path = 'sander.heat.log'

    prop = {
        "simulation_type" : "heat",
        "mdin" : { 
            'nstlim' : 2500, # Reducing the number of steps for the sake of time (5ps)
            'ntr' : 1,       # Overwritting restrain parameter
            'restraintmask' : '\"!:WAT,Cl-,Na+\"',      # Restraining solute
            'restraint_wt' : 10.0                       # With a force constant of 10 Kcal/mol*A2
        }
    }

    # Create and launch bb
    sander_mdrun(input_top_path=input_top_path, #output_ions_top_path,
                input_crd_path=input_crd_path, # output_min_rst_path,
                input_ref_path=input_crd_path, # output_min_rst_path,
                output_traj_path=output_heat_traj_path,
                output_rst_path=output_heat_rst_path,
                output_log_path=output_heat_log_path,
                properties=prop)
    
    return {
        "output_heat_traj_path": output_heat_traj_path,
        "output_heat_rst_path": output_heat_rst_path,
        "output_heat_log_path": output_heat_log_path,
    }


def equilibrate_nvt():
    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_nvt_traj_path = 'sander.nvt.netcdf'
    output_nvt_rst_path = 'sander.nvt.rst'
    output_nvt_log_path = 'sander.nvt.log'

    prop = {
        "simulation_type" : 'nvt',
        "mdin" : { 
            'nstlim' : 500, # Reducing the number of steps for the sake of time (1ps)
            'ntr' : 1,      # Overwritting restrain parameter
            'restraintmask' : '\"!:WAT,Cl-,Na+ & !@H=\"',      # Restraining solute heavy atoms
            'restraint_wt' : 5.0                               # With a force constant of 5 Kcal/mol*A2
        }
    }

    # Create and launch bb
    sander_mdrun(input_top_path=output_ions_top_path,
                input_crd_path=output_heat_rst_path,
                input_ref_path=output_heat_rst_path,
                output_traj_path=output_nvt_traj_path,
                output_rst_path=output_nvt_rst_path,
                output_log_path=output_nvt_log_path,
                properties=prop)
    
    # TODO: Clean up function with returns etc.
    

def equilibrate_npt():
    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_npt_traj_path = 'sander.npt.netcdf'
    output_npt_rst_path = 'sander.npt.rst'
    output_npt_log_path = 'sander.npt.log'

    prop = {
        "simulation_type" : 'npt',
        "mdin" : { 
            'nstlim' : 500, # Reducing the number of steps for the sake of time (1ps)
            'ntr' : 1,      # Overwritting restrain parameter
            'restraintmask' : '\"!:WAT,Cl-,Na+ & !@H=\"',      # Restraining solute heavy atoms
            'restraint_wt' : 2.5                               # With a force constant of 2.5 Kcal/mol*A2
        }
    }

    # Create and launch bb
    sander_mdrun(input_top_path=output_ions_top_path,
                input_crd_path=output_nvt_rst_path,
                input_ref_path=output_nvt_rst_path,
                output_traj_path=output_npt_traj_path,
                output_rst_path=output_npt_rst_path,
                output_log_path=output_npt_log_path,
                properties=prop)
    
    # TODO: Clean up function with returns etc.


def simulate():
    """Runs a simulation."""

    # Import module
    from biobb_amber.sander.sander_mdrun import sander_mdrun

    # Create prop dict and inputs/outputs
    output_free_traj_path = 'sander.free.netcdf'
    output_free_rst_path = 'sander.free.rst'
    output_free_log_path = 'sander.free.log'

    prop = {
        "simulation_type" : 'free',
        "mdin" : { 
            'nstlim' : 2500, # Reducing the number of steps for the sake of time (5ps)
            'ntwx' : 500  # Print coords to trajectory every 500 steps (1 ps)
        }
    }

    # Create and launch bb
    sander_mdrun(input_top_path=output_ions_top_path,
                input_crd_path=output_npt_rst_path,
                output_traj_path=output_free_traj_path,
                output_rst_path=output_free_rst_path,
                output_log_path=output_free_log_path,
                properties=prop)
    
    # TODO: Clean up function etc.
