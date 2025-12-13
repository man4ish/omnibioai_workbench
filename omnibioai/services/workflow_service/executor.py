"""
executor.py

Module for high-level orchestration of workflow execution with MySQL workflow provenance tracking.

This module defines the `WorkflowExecutor` class, which is responsible for managing
the lifecycle of workflow executions across different workflow engines
(e.g., Nextflow, Snakemake, WDL, CWL). It coordinates workflow execution, handles state
transitions, integrates with adapters, runners, and execution contexts, and persists
workflow provenance in MySQL.
"""

from typing import Iterator, Tuple
from .state import WorkflowState
from .models import WorkflowSpec
from .execution_context import ExecutionContext
from .runner import WorkflowRunner
from .adapters import (
    NextflowAdapter,
    SnakemakeAdapter,
    WDLAdapter,
    CWLAdapter 
)
from omnibioai.core.config import MYSQL_CONFIG
from .provenance import WorkflowProvenance


class WorkflowExecutor:
    """
    Orchestrates workflow execution, state transitions, and workflow provenance persistence.
    """

    ADAPTERS = {
        "nextflow": NextflowAdapter,
        "snakemake": SnakemakeAdapter,
        "wdl": WDLAdapter,
        "cwl": CWLAdapter
    }

    def execute(self, spec: WorkflowSpec) -> Iterator[Tuple[WorkflowState, int, str]]:
        """
        Execute workflow and yield state updates while persisting provenance in MySQL.

        Yields:
            (state, progress, log)
        """

        adapter_cls = self.ADAPTERS.get(spec.engine.lower())
        if not adapter_cls:
            raise ValueError(f"Unsupported engine: {spec.engine}")

        adapter = adapter_cls()
        context = ExecutionContext.create(spec.work_dir or "./work")
        runner = WorkflowRunner(adapter)

        # --------------------------
        # Initialize workflow provenance using config
        # --------------------------
        prov = WorkflowProvenance(**MYSQL_CONFIG)
        run_id = prov.create_run(spec.name, spec.engine, spec.entrypoint, spec.params)

        yield WorkflowState.RUNNING, 0, "Workflow started"

        try:
            outputs = []
            for progress, log in runner.run(spec, context):
                # Update logs in MySQL
                prov.update_progress(progress, log)
                yield WorkflowState.RUNNING, progress, log

            # Collect outputs from execution context (implement list_outputs)
            outputs = context.list_outputs()
            prov.complete_run(outputs)

            yield WorkflowState.COMPLETED, 100, "Workflow completed successfully"

        except Exception as exc:
            prov.fail_run(str(exc))
            yield WorkflowState.FAILED, 0, str(exc)
            raise

    