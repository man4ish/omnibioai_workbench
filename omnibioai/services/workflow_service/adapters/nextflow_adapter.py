"""
nextflow_adapter.py

Provides an adapter for executing workflows using the Nextflow engine.

This module defines `NextflowAdapter`, a concrete implementation of
`BaseWorkflowAdapter` that allows the OmnibioAI workflow service to
execute Nextflow pipelines via the command-line interface (CLI).

Key Features:
- CLI-based execution: Invokes Nextflow directly using subprocess,
  allowing full control over pipeline execution.
- Parameter injection: Automatically converts workflow parameters
  into CLI arguments for the Nextflow process.
- Real-time logging: Streams stdout and stderr lines as workflow
  execution progresses.
- Progress tracking: Emits incremental progress updates, useful for
  monitoring workflow status in real time.

Design Considerations:
- The adapter ensures that the working directory exists before
  execution.
- Progress reporting is heuristic and may not reflect the exact
  completion percentage of Nextflow pipelines.
- Exceptions are raised if the Nextflow process exits with a
  non-zero status.

Example Usage:
    from omnibioai.services.workflow_service.adapters.nextflow_adapter import NextflowAdapter

    adapter = NextflowAdapter(nextflow_bin="nextflow")
    for progress, log in adapter.run("main.nf", "./work", {"sample": "SAMPLE1"}):
        print(progress, log)
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
