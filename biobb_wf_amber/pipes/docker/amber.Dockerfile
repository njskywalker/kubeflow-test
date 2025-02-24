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
RUN conda init bash
COPY env/conda_env.yml /tmp/conda_env.yml
RUN conda env update -f /tmp/conda_env.yml && rm /tmp/conda_env.yml

FROM base AS final

# RUN echo "export PATH=/opt/conda/bin:$PATH" > ~/.bashrc

COPY src/stages/simulate.py /simulate.py
# in case we need to debug the container
COPY docker/mount/debug /mount/debug

# install pip and plotly (but install in source, not the conda venv)
# RUN pip install plotly

# RUN echo "conda activate amber_venv" >> ~/.bashrc
# RUN echo "source activate amber_venv" >> ~/.bashrc

# Kubeflow uses `sh` to run commands, but not through `bash`
# so `.bashrc` config is not inherited. Must create own init script
RUN echo ". /opt/conda/bin/activate" > ~/.init_sh
ENV ENV=~/.init_sh

# ENTRYPOINT [ "/bin/bash", "-l", "-c", "source activate base && conda activate base"]
CMD ["sh"]
