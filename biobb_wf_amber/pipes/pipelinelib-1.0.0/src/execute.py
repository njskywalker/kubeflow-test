import os
from kfp.client import Client

if __name__ == '__main__':
    kf_pipelines_ip = os.getenv('KUBEFLOW_PIPELINES_UI_ENDPOINT', 'http://localhost:8080')

    client = Client(host=kf_pipelines_ip)
    run = client.create_run_from_pipeline_package(
        'outputs/pipeline.yaml',  # NB Using YAML instead of passing pipe object directly
        arguments={
            'pdb_code': '1aki',
        },
    )

    print("What")
