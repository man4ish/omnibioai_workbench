from celery import Celery
from .adapters import NextflowAdapter, SnakemakeAdapter, WDLAdapter

# Celery app
app = Celery(
    "workflow_executor",
    broker="redis://localhost:6379/0",   # Adjust broker as needed
    backend="redis://localhost:6379/0"
)

# Initialize adapters
nextflow_adapter = NextflowAdapter()
snakemake_adapter = SnakemakeAdapter()
wdl_adapter = WDLAdapter()

@app.task
def run_nextflow_workflow(script_path, work_dir=None, params=None):
    return nextflow_adapter.run_workflow(script_path, work_dir, params)

@app.task
def run_snakemake_workflow(snakefile_path, work_dir=None, config_file=None):
    return snakemake_adapter.run_workflow(snakefile_path, work_dir, config_file)

@app.task
def run_wdl_workflow(wdl_file, inputs_json=None, options_json=None, work_dir=None):
    return wdl_adapter.run_workflow(wdl_file, inputs_json, options_json, work_dir)
