# start from ubuntu 20.04 image
# specify architecture of host machine
FROM --platform=linux/amd64 ubuntu:20.04 AS base

# Run as root to make sure we can install everything
# NB not for production use! :)
USER root

# install conda
RUN apt-get update && apt-get install -y wget
RUN mkdir -p ~/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    bash ~/miniconda3/miniconda.sh -b -u -p /opt/conda && \
    rm ~/miniconda3/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# create & activate conda environment
COPY env/conda_env.yml /tmp/conda_env.yml
RUN conda env create -f /tmp/conda_env.yml && rm /tmp/conda_env.yml
RUN echo "source activate biobb_wf_amber" > ~/.bashrc

# copy the pipeline code
COPY ./stages /src/stages
WORKDIR /src

FROM base AS final

# copy main code (multistage build to avoid dep installation all the time?)
COPY main.py /src/main.py

# entrypoint so I can check the conda environment
ENTRYPOINT [ "/bin/bash" ]
