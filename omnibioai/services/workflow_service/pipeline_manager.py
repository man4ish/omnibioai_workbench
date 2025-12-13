"""
pipeline_manager.py

Public API for launching and monitoring workflows.
"""

from typing import Iterator, Tuple
from .models import WorkflowSpec
from .executor import WorkflowExecutor
from .state import WorkflowState


class PipelineManager:
    """
    Entry point for workflow submission.
    """

    def __init__(self):
        self.executor = WorkflowExecutor()

    def launch(self, spec: WorkflowSpec) -> Iterator[Tuple[WorkflowState, int, str]]:
        """
        Launch a workflow and stream execution updates.
        """
        return self.executor.execute(spec)
