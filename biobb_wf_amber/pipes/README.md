## Setup

```bash
# Ensure Docker is running (macOS)
open -a Docker
# Need cluster working (on Docker) to run kf
minikube start
```

## Frustrations

1. Had to install `kfp` through `pip` instead of `conda`
as otherwise internal imports failed (`requests_toolbox` appengine)
