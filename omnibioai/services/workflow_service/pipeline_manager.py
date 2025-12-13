"""
pipeline_manager.py

Module providing the high-level public API for submitting and monitoring workflows.

This module exposes the `PipelineManager` class as the primary interface
for external users or services to launch workflows without needing to
interact directly with lower-level execution components like adapters
or runners.

Key Components:
- `PipelineManager`: Main entry point for workflow submission. Internally
  uses `WorkflowExecutor` to handle execution, state transitions, and
  progress reporting.

Responsibilities:
- Accepts declarative workflow specifications (`WorkflowSpec`) from users.
- Delegates execution to the underlying `WorkflowExecutor`.
- Streams real-time workflow updates including state, progress percentage,
  and log messages.
- Abstracts engine-specific details, allowing workflows from different engines
  (Nextflow, Snakemake, WDL) to be launched through a consistent interface.

Example Usage:
    from omnibioai.services.workflow_service.pipeline_manager import PipelineManager
    from omnibioai.services.workflow_service.models import WorkflowSpec
    from omnibioai.services.workflow_service.state import WorkflowState

    manager = PipelineManager()
    spec = WorkflowSpec(
        name="ExampleWorkflow",
        engine="nextflow",
        entrypoint="workflows/main.nf",
        params={"sample": "SAMPLE1", "threads": 4},
        work_dir="./work"
    )

    for state, progress, log in manager.launch(spec):
        print(f"{state}: {progress}% - {log}")

Purpose:
- Provide a simplified, high-level interface for workflow orchestration.
- Enable monitoring of workflow progress and status in real-time.
- Maintain separation between workflow declaration, execution, and monitoring.
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
