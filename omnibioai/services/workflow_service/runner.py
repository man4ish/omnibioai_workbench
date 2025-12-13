"""
runner.py

Core workflow execution runner module.

This module provides the `WorkflowRunner` class, which is responsible for
orchestrating the execution of workflows in a framework-agnostic manner.
It delegates the actual execution to engine-specific adapters that implement
the `BaseWorkflowAdapter` interface.

Key Components:
- `WorkflowRunner`: High-level runner that abstracts workflow execution logic
  from specific engines. Accepts a workflow adapter and executes workflows
  according to the given `WorkflowSpec` and `ExecutionContext`.

Responsibilities:
- Receive a workflow specification (`WorkflowSpec`) and execution context
  (`ExecutionContext`).
- Use the provided adapter to run the workflow script in the appropriate
  engine environment.
- Yield real-time progress updates as `(progress, log_message)` tuples.

Design Goals:
- Framework-agnostic execution: `WorkflowRunner` itself does not implement
  any engine-specific logic; it relies entirely on the adapter.
- Reusability: The runner can be reused for any workflow engine that
  implements the `BaseWorkflowAdapter` interface.
- Streaming feedback: Provides iterative updates to consumers about workflow
  progress and logs.

Example Usage:
    from omnibioai.services.workflow_service.runner import WorkflowRunner
    from omnibioai.services.workflow_service.adapters.nextflow import NextflowAdapter
    from omnibioai.services.workflow_service.models import WorkflowSpec
    from omnibioai.services.workflow_service.execution_context import ExecutionContext

    adapter = NextflowAdapter()
    runner = WorkflowRunner(adapter)
    spec = WorkflowSpec(
        name="ExampleWorkflow",
        engine="nextflow",
        entrypoint="workflows/main.nf",
        params={"sample": "SAMPLE1", "threads": 4},
        work_dir="./work"
    )
    context = ExecutionContext.create(spec.work_dir)

    for progress, log in runner.run(spec, context):
        print(f"{progress}% - {log}")
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
