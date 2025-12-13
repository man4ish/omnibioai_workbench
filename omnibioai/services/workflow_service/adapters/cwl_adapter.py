"""
cwl_adapter.py

Provides an adapter for executing workflows using the Common Workflow Language (CWL).

This module defines `CWLAdapter`, a concrete implementation of
`BaseWorkflowAdapter` that allows the OmnibioAI workflow service to
execute CWL workflows via the `cwltool` command-line interface.

Key Features:
- CLI-based execution: Invokes cwltool directly using subprocess.
- Parameter injection: Converts workflow inputs into CLI arguments.
- Real-time logging: Streams stdout and stderr lines as workflow
  execution progresses.
- Progress tracking: Emits incremental progress updates, useful for
  monitoring workflow status in real time.

Design Considerations:
- The adapter ensures that the working directory exists before execution.
- Progress reporting is heuristic and may not reflect exact workflow
  completion percentage.
- Exceptions are raised if the CWL process exits with a non-zero status.

Example Usage:
    from omnibioai.services.workflow_service.adapters.cwl_adapter import CWLAdapter

    adapter = CWLAdapter(cwltool_bin="cwltool")
    for progress, log in adapter.run("workflow.cwl", "./work", {"input_file": "sample.fastq"}):
        print(progress, log)
"""

import subprocess
import os
from typing import Iterator, Tuple
from .base import BaseWorkflowAdapter


class CWLAdapter(BaseWorkflowAdapter):
    """
    Executes CWL workflows via cwltool CLI.
    """

    def __init__(self, cwltool_bin: str = "cwltool"):
        self.cwltool_bin = cwltool_bin

    def run(
        self,
        entrypoint: str,
        work_dir: str,
        params: dict
    ) -> Iterator[Tuple[int, str]]:

        os.makedirs(work_dir, exist_ok=True)

        cmd = [self.cwltool_bin, entrypoint]
        for k, v in params.items():
            cmd.append(f"--{k}")
            cmd.append(str(v))

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
            raise RuntimeError("CWL execution failed")

        yield 100, "CWL execution completed"

