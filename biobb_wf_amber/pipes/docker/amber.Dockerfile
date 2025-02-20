# FROM --platform=linux/amd64 condaforge/mambaforge:latest AS base
FROM --platform=linux/amd64 quay.io/biocontainers/biobb_amber:5.0.4--pyhdfd78af_0 AS base
SHELL ["/bin/bash", "--login", "-c"]

# grab and install conda
RUN mkdir -p ~/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    bash ~/miniconda3/miniconda.sh -b -u -p /opt/conda && \
    rm ~/miniconda3/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# create & activate conda environment
COPY env/conda_env.yml /tmp/conda_env.yml
RUN conda env create -f /tmp/conda_env.yml && rm /tmp/conda_env.yml
RUN echo "conda activate biobb_wf_amber" > ~/.bashrc

FROM base AS final

RUN conda init bash
ENV PATH=/opt/conda/bin:$PATH

# ENTRYPOINT [ "/bin/bash", "--login", "-c", "conda activate biobb_wf_amber" ]
ENTRYPOINT [ "/bin/bash", "--login", "-c", "conda init bash && . ~/.bashrc && conda activate biobb_wf_amber"]
# ENTRYPOINT ["conda", "run", "-n", "biobb_wf_amber", "python3"]
