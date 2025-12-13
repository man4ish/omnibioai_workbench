"""
nextflow_adapter.py

Adapter for executing Nextflow workflows.
"""

import subprocess
import os
from typing import Iterator, Tuple
from .base import BaseWorkflowAdapter


class NextflowAdapter(BaseWorkflowAdapter):
    """
    Executes Nextflow pipelines via CLI.
    """

    def __init__(self, nextflow_bin: str = "nextflow"):
        self.nextflow_bin = nextflow_bin

    def run(
        self,
        entrypoint: str,
        work_dir: str,
        params: dict
    ) -> Iterator[Tuple[int, str]]:

        os.makedirs(work_dir, exist_ok=True)

        cmd = [self.nextflow_bin, "run", entrypoint]
        for k, v in params.items():
            cmd.append(f"--{k}={v}")

        process = subprocess.Popen(
            cmd,
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        progress = 0
        for line in process.stdout:
            progress = min(progress + 1, 100)
            yield progress, line.strip()

        if process.wait() != 0:
            raise RuntimeError("Nextflow execution failed")

        yield 100, "Nextflow execution completed"
