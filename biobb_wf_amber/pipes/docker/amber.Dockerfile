# FROM --platform=linux/amd64 condaforge/mambaforge:latest AS base
FROM --platform=linux/amd64 quay.io/biocontainers/biobb_amber:5.0.4--pyhdfd78af_0 AS base
SHELL ["/bin/bash", "--login", "-c"]

# grab and install conda
RUN mkdir -p ~/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    bash ~/miniconda3/miniconda.sh -b -u -p /opt/conda && \
    rm ~/miniconda3/miniconda.sh
RUN echo "export PATH=/opt/conda/bin:$PATH" > ~/.bashrc

# create & activate conda environment
COPY env/conda_env.yml /tmp/conda_env.yml
RUN conda env update -f /tmp/conda_env.yml && rm /tmp/conda_env.yml

FROM base AS final

RUN echo "export PATH=/opt/conda/bin:$PATH" > ~/.bashrc
RUN conda init bash

COPY src/stages/simulate.py /simulate.py

# in case we need to debug the container
COPY docker/mount/debug /mount/debug

# ENTRYPOINT [ "/bin/bash", "--login", "-c", "conda init bash && . ~/.bashrc && source activate biobb_wf_amber && conda info --envs"]
ENTRYPOINT [ "/bin/bash", "--login", "-c"]
