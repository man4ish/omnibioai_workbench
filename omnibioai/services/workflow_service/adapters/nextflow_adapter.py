import subprocess
import os

class NextflowAdapter:
    """
    Adapter to run Nextflow pipelines
    """

    def __init__(self, nextflow_bin="nextflow"):
        self.nextflow_bin = nextflow_bin

    def run_pipeline(self, script, params, work_dir=None):
        """
        Run Nextflow pipeline with optional params dictionary
        """
        cmd = [self.nextflow_bin, "run", script]

        if params:
            for k, v in params.items():
                cmd.append(f"--{k}={v}")

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            cmd.extend(["-w", work_dir])

        process = subprocess.run(cmd, capture_output=True, text=True)
        if process.returncode != 0:
            raise RuntimeError(f"Nextflow pipeline failed: {process.stderr}")
        return process.stdout
