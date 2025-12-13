"""
base.py

Defines the abstract base class for workflow engine adapters in the
OmnibioAI workflow service.

This module provides the `BaseWorkflowAdapter` class, which specifies
the interface that all concrete workflow engine adapters must implement.
Adapters act as a bridge between the high-level workflow execution
logic and specific workflow engines such as Nextflow, Snakemake, or WDL.

Key Components:
- `BaseWorkflowAdapter`: Abstract base class that defines the required
  methods for all workflow engine adapters. Ensures a consistent
  interface for executing workflows regardless of the underlying engine.

Design Goals:
- Engine-agnostic execution: Decouples workflow orchestration from
  specific workflow engines.
- Extensibility: Supports addition of new workflow engines by creating
  adapters that implement the base interface.
- Streaming execution: Adapters yield incremental progress updates
  and log messages to allow real-time monitoring.

Example Usage:
    from omnibioai.services.workflow_service.adapters.base import BaseWorkflowAdapter

    class DummyAdapter(BaseWorkflowAdapter):
        def run(self, entrypoint: str, work_dir: str, params: dict):
            for i in range(0, 101, 10):
                yield i, f"Progress: {i}%"
"""


from abc import ABC, abstractmethod
from typing import Iterator, Tuple


class BaseWorkflowAdapter(ABC):
    """
    Base adapter contract for workflow engines.
    """

    @abstractmethod
    def run(
        self,
        entrypoint: str,
        work_dir: str,
        params: dict
    ) -> Iterator[Tuple[int, str]]:
        """
        Execute the workflow.

        Yields:
            (progress_percent, log_line)

        Notes:
            Progress may be heuristic.
        """
        pass
