"""Functions for processing MD outputs for analysis."""

def process_outs():
    # Import module
    from biobb_amber.process.process_mdout import process_mdout

    # Create prop dict and inputs/outputs
    output_dat_npt_path = 'sander.md.npt.dat'

    prop = {
        "terms" : ['PRES','DENSITY']
    }

    # Create and launch bb
    process_mdout(input_log_path=output_npt_log_path,
                output_dat_path=output_dat_npt_path,
                properties=prop)