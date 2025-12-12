import subprocess
import os

class SnakemakeAdapter:
    """
    Adapter to run Snakemake workflows
    """

    def __init__(self, snakemake_bin="snakemake"):
        self.snakemake_bin = snakemake_bin

    def run_workflow(self, workflow_file, config=None, work_dir=None):
        """
        Run Snakemake workflow
        """
        cmd = [self.snakemake_bin, "--snakefile", workflow_file]

        if config:
            for k, v in config.items():
                cmd.extend(["--config", f"{k}={v}"])

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            cmd.extend(["--directory", work_dir])

        process = subprocess.run(cmd, capture_output=True, text=True)
        if process.returncode != 0:
            raise RuntimeError(f"Snakemake workflow failed: {process.stderr}")
        return process.stdout
