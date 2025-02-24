.PHONY: start-kf package build compile compile-only run pf-kf pf-minio

SHELL := /bin/bash
export PIPELINE_VERSION:=2.3.0
export KUBEFLOW_PIPELINES_UI_ENDPOINT:=http://localhost:8080
export KFP_PLATFORM:=platform-agnostic-emissary


setup-minikube:
	# Start minikube
	minikube start --cpus=4 --memory=8192
	minikube addons enable registry

	# Install kubeflow pipelines
	kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/cluster-scoped-resources?ref=$(PIPELINE_VERSION)"
	kubectl wait --for condition=established --timeout=60s crd/applications.app.k8s.io
	kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/env/$(KFP_PLATFORM)?ref=$(PIPELINE_VERSION)"

build:
	eval $$(minikube docker-env -u) && docker build -f docker/base.Dockerfile -t pipeline:latest .

compile-only:
	eval $$(minikube docker-env -u) && \
	docker run \
		--rm --network host \
		-v ./outputs:/src/outputs \
		pipeline:latest

compile: build compile-only

# TODO Do not use uv
run:
	uv run python src/execute.py

pf-kf:
	kubectl port-forward -n kubeflow svc/ml-pipeline-ui 8080:80

pf-registry:
	kubectl port-forward --namespace kube-system service/registry 5000:80

pf-minio:
	kubectl port-forward -n kubeflow svc/minio-service 9000:9000


# AMBER image (nebjovanovic/amber_bio:latest)
build-amber:
	eval $$(minikube docker-env) && docker build -f docker/amber.Dockerfile -t amber_bio:latest .

debug-mount-amber:
	docker run -it \
		-v /Users/neb.jovanovic/Documents/projects/kubeflow-test/biobb_wf_amber/pipes/docker/mount:/mount \
		--entrypoint /bin/bash amber_bio:latest
