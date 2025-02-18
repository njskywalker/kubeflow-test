## Setup

Testing commit.

```bash
# Install Docker, minikube, kubectx
sudo snap install kubectx --classic

# Ensure Docker is running (macOS)
open -a Docker

# Need cluster working (on Docker) to run kf
minikube start

# Deploy kfp resources
export PIPELINE_VERSION=2.3.0
kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/cluster-scoped-resources?ref=$PIPELINE_VERSION"
kubectl wait --for condition=established --timeout=60s crd/applications.app.k8s.io
kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/env/dev?ref=$PIPELINE_VERSION"

```

## Experience

### MinIO
I chose to use this as a local storage server in between pipeline steps. It
comes with Kubeflow and is relatively straightforward to set up and use for
local demonstration.

You can access the MinIO UI by port forwarding to the service and logging in
with the `minio/minio123` default credentials:

```bash
make port-forward-minio
```

Optionally, install the CLI to access MinIO directly from terminal, instead of
port forwarding  (Brew or Snap for macOS / Linux respectively). Then, configure
to talk to KFlow's instance:

```bash
mc config host add minio http://localhost:9000 minio minio123
```

## Frustrations

1. Had to install `kfp` through `pip` instead of `conda`
as otherwise internal imports failed (`requests_toolbox` appengine)
2. O God. Never use WSL. Just invest in a Mac
3. 

## Choices

### Infrastructure
Options considered:
1. `minikube`
2. GKE

Using GKE would require setting up the infra (e.g. with Terraform or similar,
or ClickOps) 

### Storage
Options considered:
1. MinIO
2. Volumes
3. Direct input / output passing

