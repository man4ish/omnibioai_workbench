"""
execution_context.py

Module for managing runtime metadata of workflow executions.

This module provides the `ExecutionContext` data class, which encapsulates
essential metadata for a single workflow run. It is designed to be used
by workflow adapters and pipeline managers to track, organize, and manage
workflow executions consistently across different workflow engines
(e.g., Nextflow, Snakemake, WDL, CWL).

The ExecutionContext includes:
- `run_id`: A unique identifier for each workflow execution.
- `start_time`: The UTC timestamp marking the start of execution.
- `work_dir`: The working directory where workflow outputs and logs are stored.

Key Features:
- Provides a standardized representation of workflow execution metadata.
- Includes a static factory method `create()` to generate new execution contexts
  with automatically assigned unique IDs and start timestamps.
- Facilitates integration with logging, monitoring, and workflow state tracking
  components of the OmniBioAI workflow service.
- Adds `list_outputs()` method to collect workflow output files.
"""

from dataclasses import dataclass
from datetime import datetime
import uuid
import os


@dataclass
class ExecutionContext:
    """
    Runtime execution context.

    Attributes:
        run_id: Unique execution identifier
        start_time: Execution start timestamp
        work_dir: Execution working directory
    """
    run_id: str
    start_time: datetime
    work_dir: str

    @staticmethod
    def create(work_dir: str) -> "ExecutionContext":
        os.makedirs(work_dir, exist_ok=True)
        return ExecutionContext(
            run_id=str(uuid.uuid4()),
            start_time=datetime.utcnow(),
            work_dir=work_dir
        )

    def list_outputs(self):
        """
        Return a list of all files in the workflow working directory,
        including subdirectories. Can be used to persist workflow outputs.
        """
        outputs = []
        for root, dirs, files in os.walk(self.work_dir):
            for file in files:
                outputs.append(os.path.join(root, file))
        return outputs
