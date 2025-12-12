from .executor import run_nextflow_workflow, run_snakemake_workflow, run_wdl_workflow
from celery.result import AsyncResult

class PipelineManager:
    """
    Central manager to schedule, track, and launch pipelines
    """

    def __init__(self):
        self.active_pipelines = {}

    def launch_nextflow(self, script_path, user_id, work_dir=None, params=None):
        task = run_nextflow_workflow.delay(script_path, user_id, work_dir, params)
        self.active_pipelines[task.id] = {"type": "nextflow", "status": "queued"}
        return task.id

    def launch_snakemake(self, snakefile_path, user_id, work_dir=None, config_file=None):
        task = run_snakemake_workflow.delay(snakefile_path, user_id, work_dir, config_file)
        self.active_pipelines[task.id] = {"type": "snakemake", "status": "queued"}
        return task.id

    def launch_wdl(self, wdl_file, user_id, inputs_json=None, options_json=None, work_dir=None):
        task = run_wdl_workflow.delay(wdl_file, user_id, inputs_json, options_json, work_dir)
        self.active_pipelines[task.id] = {"type": "wdl", "status": "queued"}
        return task.id

    def get_status(self, task_id):
        result = AsyncResult(task_id)
        return {"status": result.status, "result": result.result}
