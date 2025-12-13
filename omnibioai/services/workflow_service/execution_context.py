"""
execution_context.py

Captures runtime metadata for a workflow execution.
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
