FROM --platform=linux/amd64 python:3.10-slim AS base
# FROM --platform=linux/amd64 condaforge/mambaforge:latest AS base
SHELL ["/bin/bash", "-c"]

# install 
RUN apt-get update
# RUN mkdir -p ~/miniconda3 && \
#     wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
#     bash ~/miniconda3/miniconda.sh -b -u -p /opt/conda && \
#     rm ~/miniconda3/miniconda.sh
# ENV PATH=/opt/conda/bin:$PATH

# create & activate conda environment
# COPY env/conda_env.yml /tmp/conda_env.yml
# RUN conda env create -f /tmp/conda_env.yml && rm /tmp/conda_env.yml
# RUN echo "source activate biobb_wf_amber" > ~/.bashrc

RUN pip install kfp kfp-kubernetes

FROM base AS final
SHELL ["/bin/bash", "-c"]

# copy the pipeline code
COPY ./src /src
WORKDIR /src

# copy custom package
# need classes so Kubeflow can find them
# COPY dist/ /src/dist/
# RUN source activate biobb_wf_amber && pip install /src/dist/pipelinelib-1.0.0-py3-none-any.whl

COPY ./src/pipeline.py pipeline.py
ENTRYPOINT ["python", "pipeline.py"]
