"""
wdl_adapter.py

Adapter for executing WDL workflows using Cromwell.
"""

import subprocess
import os
from typing import Iterator, Tuple
from .base import BaseWorkflowAdapter


class WDLAdapter(BaseWorkflowAdapter):
    """
    Executes WDL workflows via Cromwell.
    """

    def __init__(self, cromwell_jar: str = "cromwell.jar", java_bin: str = "java"):
        self.cromwell_jar = cromwell_jar
        self.java_bin = java_bin

    def run(
        self,
        entrypoint: str,
        work_dir: str,
        params: dict
    ) -> Iterator[Tuple[int, str]]:

        os.makedirs(work_dir, exist_ok=True)

        cmd = [
            self.java_bin, "-jar", self.cromwell_jar,
            "run", entrypoint
        ]

        if "inputs" in params:
            cmd.extend(["-i", params["inputs"]])
        if "options" in params:
            cmd.extend(["-o", params["options"]])

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
            raise RuntimeError("WDL execution failed")

        yield 100, "WDL execution completed"
