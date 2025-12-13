"""
runner.py

Core workflow execution runner.
Framework-agnostic and reusable.
"""

from typing import Iterator, Tuple
from .models import WorkflowSpec
from .execution_context import ExecutionContext
from .adapters.base import BaseWorkflowAdapter


class WorkflowRunner:
    """
    Executes a workflow using a specific adapter.
    """

    def __init__(self, adapter: BaseWorkflowAdapter):
        self.adapter = adapter

    def run(
        self,
        spec: WorkflowSpec,
        context: ExecutionContext
    ) -> Iterator[Tuple[int, str]]:

        yield from self.adapter.run(
            entrypoint=spec.entrypoint,
            work_dir=context.work_dir,
            params=spec.params
        )
