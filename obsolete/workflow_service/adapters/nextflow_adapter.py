import subprocess
import os

class NextflowAdapter:
    """
    Adapter to run Nextflow pipelines with progress reporting
    """

    def __init__(self, nextflow_bin="nextflow"):
        self.nextflow_bin = nextflow_bin

    def run_workflow(self, script, work_dir=None, params=None, progress=False):
        """
        Run Nextflow pipeline with optional params dictionary.
        Yields (progress, log) if progress=True
        """
        cmd = [self.nextflow_bin, "run", script]

        if params:
            for k, v in params.items():
                cmd.append(f"--{k}={v}")

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            cmd.extend(["-w", work_dir])

        # Run Nextflow as subprocess
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        if progress:
            step = 0
            for line in process.stdout:
                line = line.strip()
                step += 1
                yield min(step, 100), line  # simplistic progress
        else:
            output, _ = process.communicate()
            if process.returncode != 0:
                raise RuntimeError(f"Nextflow pipeline failed: {output}")
            yield 100, output
