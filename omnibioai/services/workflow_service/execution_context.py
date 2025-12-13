"""
execution_context.py

Module for managing runtime metadata of workflow executions.

This module provides the `ExecutionContext` data class, which encapsulates
essential metadata for a single workflow run. It is designed to be used
by workflow adapters and pipeline managers to track, organize, and manage
workflow executions consistently across different workflow engines
(e.g., Nextflow, Snakemake, WDL).

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

Example Usage:
    from omnibioai.services.workflow_service.execution_context import ExecutionContext

    # Create a new execution context for a workflow
    ctx = ExecutionContext.create(work_dir="/path/to/workdir")

    print(ctx.run_id)       # Unique execution ID
    print(ctx.start_time)   # Execution start timestamp
    print(ctx.work_dir)     # Working directory path
"""


from dataclasses import dataclass
from datetime import datetime
import uuid


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
        return ExecutionContext(
            run_id=str(uuid.uuid4()),
            start_time=datetime.utcnow(),
            work_dir=work_dir
        )
