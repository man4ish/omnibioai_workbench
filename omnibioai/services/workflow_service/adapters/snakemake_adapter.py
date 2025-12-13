"""
snakemake_adapter.py

Provides an adapter for executing workflows using the Snakemake engine.

This module defines `SnakemakeAdapter`, a concrete implementation that
allows the OmnibioAI workflow service to execute Snakemake pipelines
via the command-line interface (CLI).

Key Features:
- CLI-based execution: Invokes Snakemake directly using subprocess,
  supporting full pipeline control.
- Configurable execution: Supports `cores`, `config` dictionaries,
  working directory, and additional parameters.
- Real-time logging: Streams stdout and stderr lines as the workflow
  progresses when `progress=True`.
- Progress reporting: Emits incremental step-based progress updates
  to help monitor workflow execution.

Design Considerations:
- Ensures that the working directory exists before execution.
- Progress is heuristic and may not reflect exact workflow completion.
- Raises exceptions if the Snakemake process exits with a non-zero
  status.

Example Usage:
    from omnibioai.services.workflow_service.adapters.snakemake_adapter import SnakemakeAdapter

    adapter = SnakemakeAdapter(snakemake_bin="snakemake")
    for progress, log in adapter.run(
        entrypoint="Snakefile",
        work_dir="./work",
        config={"sample": "SAMPLE1"},
        progress=True,
        params={"cores": 4}
    ):
        print(progress, log)
"""

import subprocess
import os

class SnakemakeAdapter:
    """
    Adapter to run Snakemake workflows with progress reporting
    """

    def __init__(self, snakemake_bin="snakemake"):
        self.snakemake_bin = snakemake_bin

    def run(self, entrypoint, work_dir=None, config=None, progress=False, params=None):
        """
        Run Snakemake workflow. Yields (progress, log) if progress=True.

        Parameters
        ----------
        entrypoint : str
            Path to Snakefile
        work_dir : str
            Working directory
        config : dict
            Config dictionary for Snakemake
        progress : bool
            Yield step progress
        params : dict
            Extra params like {"cores": 1}
        """
        cmd = [self.snakemake_bin, "--snakefile", entrypoint, "--printshellcmds", "--rerun-incomplete"]

        # Set cores
        cores = 1
        if params and "cores" in params:
            cores = params["cores"]
        cmd.extend(["--cores", str(cores)])

        # Config
        if config:
            for k, v in config.items():
                cmd.extend(["--config", f"{k}={v}"])

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            cmd.extend(["--directory", work_dir])

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        if progress:
            step = 0
            for line in process.stdout:
                line = line.strip()
                step += 1
                yield min(step, 100), line
        else:
            output, _ = process.communicate()
            if process.returncode != 0:
                raise RuntimeError(f"Snakemake workflow failed: {output}")
            yield 100, output


