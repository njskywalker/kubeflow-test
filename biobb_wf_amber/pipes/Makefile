.PHONY: start-kf package build compile compile-only run pf-kf pf-minio

SHELL := /bin/bash
export PIPELINE_VERSION:=2.3.0
export KUBEFLOW_PIPELINES_UI_ENDPOINT:=http://localhost:8080
export KFP_PLATFORM:=platform-agnostic-emissary


setup-kf:
	minikube addons enable registry  # For Dockerimages 
	kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/cluster-scoped-resources?ref=$(PIPELINE_VERSION)"
	kubectl wait --for condition=established --timeout=60s crd/applications.app.k8s.io
	kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/env/$(KFP_PLATFORM)?ref=$(PIPELINE_VERSION)"

package:
	pip install --upgrade setuptools wheel
	python setup.py sdist bdist_wheel
	pip install dist/pipelinelib-1.0.0.tar.gz

build:
	docker build -f docker/base.Dockerfile -t localhost:5000/pipeline:latest .
	docker push localhost:5000/pipeline:latest

compile-only:
	docker run -it \
		--rm --network host \
		-v ./outputs:/src/outputs \
		pipeline:latest

compile: build compile-only

run:
	uv run python src/execute.py

pf-kf:
	kubectl port-forward -n kubeflow svc/ml-pipeline-ui 8080:80

pf-minio:
	kubectl port-forward -n kubeflow svc/minio-service 9000:9000
