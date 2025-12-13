"""
wdl_adapter.py

Provides an adapter for executing workflows written in the Workflow Description
Language (WDL) using the Cromwell execution engine.

This module defines `WDLAdapter`, a concrete implementation that allows
OmnibioAI to run WDL pipelines via the Cromwell CLI. It provides optional
support for input JSON files and Cromwell runtime options, as well as
streaming execution logs and progress updates.

Key Features:
- CLI-based execution: Runs WDL workflows by invoking Cromwell via Java.
- Configurable execution: Supports `inputs` and `options` via parameters.
- Real-time logging: Streams stdout and stderr as the workflow runs.
- Progress reporting: Emits heuristic progress updates to monitor execution.
- Automatic working directory setup.

Design Considerations:
- Ensures the working directory exists before workflow execution.
- Progress is step-based and heuristic; it may not reflect precise workflow completion.
- Raises exceptions if the Cromwell process exits with a non-zero status.

Example Usage:
    from omnibioai.services.workflow_service.adapters.wdl_adapter import WDLAdapter

    adapter = WDLAdapter(cromwell_jar="cromwell.jar", java_bin="java")
    for progress, log in adapter.run(
        entrypoint="workflow.wdl",
        work_dir="./work",
        params={
            "inputs": "inputs.json",
            "options": "options.json"
        }
    ):
        print(progress, log)
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
