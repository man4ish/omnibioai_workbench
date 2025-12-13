"""
executor.py

High-level workflow executor that manages lifecycle state.
"""

from typing import Iterator, Tuple
from .state import WorkflowState
from .models import WorkflowSpec
from .execution_context import ExecutionContext
from .runner import WorkflowRunner
from .adapters import (
    NextflowAdapter,
    SnakemakeAdapter,
    WDLAdapter
)


class WorkflowExecutor:
    """
    Orchestrates workflow execution and state transitions.
    """

    ADAPTERS = {
        "nextflow": NextflowAdapter,
        "snakemake": SnakemakeAdapter,
        "wdl": WDLAdapter
    }

    def execute(self, spec: WorkflowSpec) -> Iterator[Tuple[WorkflowState, int, str]]:
        """
        Execute workflow and yield state updates.

        Yields:
            (state, progress, log)
        """

        adapter_cls = self.ADAPTERS.get(spec.engine.lower())
        if not adapter_cls:
            raise ValueError(f"Unsupported engine: {spec.engine}")

        adapter = adapter_cls()
        context = ExecutionContext.create(spec.work_dir or "./work")

        runner = WorkflowRunner(adapter)

        yield WorkflowState.RUNNING, 0, "Workflow started"

        try:
            for progress, log in runner.run(spec, context):
                yield WorkflowState.RUNNING, progress, log

            yield WorkflowState.COMPLETED, 100, "Workflow completed successfully"

        except Exception as exc:
            yield WorkflowState.FAILED, 0, str(exc)
            raise
