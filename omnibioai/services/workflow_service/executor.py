"""
executor.py

Module for high-level orchestration of workflow execution.

This module defines the `WorkflowExecutor` class, which is responsible for managing
the lifecycle of workflow executions across different workflow engines
(e.g., Nextflow, Snakemake, WDL). It coordinates workflow execution, handles state
transitions, and integrates with adapters, runners, and execution contexts.

Key Components:
- `WorkflowExecutor`: The main class that executes workflows and yields real-time
  state updates, progress, and log messages.
- `ADAPTERS`: A mapping of supported workflow engines to their respective adapter
  classes. Adapters provide engine-specific execution capabilities.
- `execute()`: Method that orchestrates the execution of a workflow based on a
  `WorkflowSpec`, yielding tuples of `(WorkflowState, progress, log)`.

Features:
- Supports multiple workflow engines through adapter classes.
- Generates an `ExecutionContext` for each workflow run to track metadata
  such as run ID, start time, and working directory.
- Provides real-time feedback via a generator yielding workflow state,
  progress percentage, and log messages.
- Handles exceptions gracefully, marking the workflow as `FAILED` if an error occurs.

Example Usage:
    from omnibioai.services.workflow_service.executor import WorkflowExecutor
    from omnibioai.services.workflow_service.models import WorkflowSpec

    spec = WorkflowSpec(engine="nextflow", work_dir="./demo")
    executor = WorkflowExecutor()

    for state, progress, log in executor.execute(spec):
        print(f"{state.name}: {progress}% - {log}")

    # Output will stream workflow progress and final state
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
