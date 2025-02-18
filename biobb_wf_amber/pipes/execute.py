import os
from kfp.client import Client

minikube_ip = os.getenv('MINIKUBE_IP', 'localhost:8080')

client = Client(host=minikube_ip)
run = client.create_run_from_pipeline_package(
    'outputs/pipeline.yaml',  # NB Using YAML instead of passing pipe object directly
    arguments={
        'pdb_code': '1aki',
    },
)
