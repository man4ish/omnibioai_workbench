import time
from omnibioai.core.services.workflow_service import PipelineManager

def main():
    pm = PipelineManager()

    # ----------------------------
    # Launch Nextflow workflow
    # ----------------------------
    nextflow_script = "workflows/example.nf"
    nextflow_params = {"sample": "SAMPLE1", "threads": 4}
    nextflow_task_id = pm.launch_nextflow(nextflow_script, work_dir="workflows/nf_workdir", params=nextflow_params)
    print(f"Nextflow workflow queued with task ID: {nextflow_task_id}")

    # ----------------------------
    # Launch Snakemake workflow
    # ----------------------------
    snakefile = "workflows/Snakefile"
    snakemake_task_id = pm.launch_snakemake(snakefile, work_dir="workflows/snake_workdir", config_file="workflows/config.yaml")
    print(f"Snakemake workflow queued with task ID: {snakemake_task_id}")

    # ----------------------------
    # Launch WDL workflow
    # ----------------------------
    wdl_file = "workflows/example.wdl"
    inputs_json = "workflows/inputs.json"
    options_json = "workflows/options.json"
    wdl_task_id = pm.launch_wdl(wdl_file, inputs_json, options_json, work_dir="workflows/wdl_workdir")
    print(f"WDL workflow queued with task ID: {wdl_task_id}")

    # ----------------------------
    # Poll pipeline statuses
    # ----------------------------
    task_ids = [nextflow_task_id, snakemake_task_id, wdl_task_id]
    while task_ids:
        for task_id in task_ids[:]:
            status_info = pm.get_status(task_id)
            print(f"Task {task_id} | Type: {pm.active_pipelines[task_id]['type']} | Status: {status_info['status']}")
            if status_info['status'] in ["SUCCESS", "FAILURE", "REVOKED"]:
                print(f"Result for task {task_id}: {status_info['result']}")
                task_ids.remove(task_id)
        time.sleep(5)

if __name__ == "__main__":
    main()
