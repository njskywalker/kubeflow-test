FROM amber_bio:latest AS base
SHELL ["/bin/bash", "--login", "-c"]

COPY src/analysis/simulation.py /analysis/simulation.py
RUN pip install plotly

ENTRYPOINT [ "/bin/bash", "--login", "-c"]
