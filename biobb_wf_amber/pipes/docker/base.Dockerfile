FROM python:3.10-slim AS base
SHELL ["/bin/bash", "-c"]

# update & install deps
RUN apt-get update
RUN pip install kfp kfp-kubernetes

FROM base AS final
SHELL ["/bin/bash", "-c"]

# copy the pipeline code
COPY ./src /src
WORKDIR /src

ENTRYPOINT ["python", "pipeline.py"]
