""""""

from kfp import dsl
from typing import Literal

@dsl.component(
    base_image="amber_bio:latest",
    # need install plotly like this
    # again, conda env doesn't load in image
    packages_to_install=["plotly"],
)
def plot_simulation_scatter(
    properties: dict,
    process_type: str,
    plot_title: str,
    x_axis_title: str,
    y_axis_title: str,
    plot_height: int,
    # kubeflow 
    input_path: dsl.InputPath('Directory'),
    input_log_filename: str,
    output_path: dsl.OutputPath('Directory'),
    output_dat_filename: str,
):
    """Process simulation log files to generate and save 
    a scatter plot."""

    # Import modules
    import os
    import plotly.graph_objs as go
    if process_type == "minout":
        from biobb_amber.process.process_minout import process_minout
        process_func = process_minout
    else:
        from biobb_amber.process.process_mdout import process_mdout
        process_func = process_mdout

    # Process log files, creating stats
    os.makedirs(output_path, exist_ok=True)
    input_log_path = input_path + "/" + input_log_filename
    output_dat_path = output_path + "/" + output_dat_filename


    # Create and launch bb
    process_func(input_log_path=input_log_path,
                output_dat_path=output_dat_path,
                properties=properties)  

    with open(output_dat_path, 'r') as energy_file:
        x, y = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in energy_file
            if not line.startswith(("#", "@"))
            if float(line.split()[1]) < 1000
        ])

    # Create a scatter plot
    fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

    # Update layout
    fig.update_layout(title=plot_title,
                    xaxis_title=x_axis_title,
                    yaxis_title=y_axis_title,
                    height=plot_height)

    # Save the plot
    fig.write_html(output_path + "/simulation_scatter_plot.html")
