.PHONY: start-kf

SHELL := /bin/bash
export PIPELINE_VERSION:=2.3.0
export MINIKUBE_IP:=http://$(shell minikube ip):8080

build:
	docker build -f docker/base.Dockerfile -t pipeline:latest .

compile:
	docker run -it \
		--rm --network host \
		-v ./outputs:/src/outputs \
		pipeline:latest

start-kf:
	kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/cluster-scoped-resources?ref=$(PIPELINE_VERSION)"
	kubectl wait --for condition=established --timeout=60s crd/applications.app.k8s.io
	kubectl apply -k "github.com/kubeflow/pipelines/manifests/kustomize/env/dev?ref=$(PIPELINE_VERSION)"

port-forward-kf:
	kubectl port-forward -n kubeflow svc/ml-pipeline-ui 8080:80

port-forward-minio:
	kubectl port-forward -n kubeflow svc/minio-service 9000:9000
