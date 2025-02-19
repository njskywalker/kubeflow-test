# FROM --platform=linux/amd64 ubuntu:20.04 AS base
FROM --platform=linux/amd64 condaforge/mambaforge:latest AS base
SHELL ["/bin/bash", "-c"]

# # Run as root to make sure we can install everything
# # NB not for production use! :)
# USER root

# install conda
RUN apt-get update
# RUN mkdir -p ~/miniconda3 && \
#     wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
#     bash ~/miniconda3/miniconda.sh -b -u -p /opt/conda && \
#     rm ~/miniconda3/miniconda.sh
# ENV PATH=/opt/conda/bin:$PATH

# create & activate conda environment
COPY env/conda_env.yml /tmp/conda_env.yml
RUN conda env create -f /tmp/conda_env.yml && rm /tmp/conda_env.yml
RUN echo "source activate biobb_wf_amber" > ~/.bashrc

# copy the pipeline code
COPY ./src /src
WORKDIR /src

FROM base AS final

# copy custom package
# need classes so Kubeflow can find them
COPY dist/ /src/dist/
SHELL ["/bin/bash", "-c"]
RUN source activate biobb_wf_amber && pip install /src/dist/pipelinelib-1.0.0-py3-none-any.whl

COPY ./src/compile.py compile.py
ENTRYPOINT [ "/bin/bash", "-c", "source activate biobb_wf_amber && exec python compile.py" ]
