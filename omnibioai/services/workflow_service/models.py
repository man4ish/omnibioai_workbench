"""
models.py

Module defining core data models for workflow execution.

This module contains framework-agnostic data classes that describe workflows
declaratively. These models do not contain execution logic themselves, but
serve as structured representations of workflows that can be consumed by
workflow executors, adapters, and runners.

Key Classes:
- `WorkflowSpec`: Encapsulates the full specification of a workflow, including
  its name, engine, entrypoint script, execution parameters, and optional
  working directory. This class is used as the main input to the
  `WorkflowExecutor`.

Purpose:
- Provides a consistent and type-safe way to define workflows across
  different workflow engines (Nextflow, Snakemake, WDL).
- Enables separation of workflow **declaration** from **execution logic**.
- Facilitates code readability, validation, and maintainability.

Example Usage:
    from omnibioai.services.workflow_service.models import WorkflowSpec

    spec = WorkflowSpec(
        name="ExampleWorkflow",
        engine="nextflow",
        entrypoint="workflows/main.nf",
        params={"sample": "SAMPLE1", "threads": 4},
        work_dir="./work"
    )

    # `spec` can now be passed to a WorkflowExecutor for execution
"""


from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class WorkflowSpec:
    """
    Declarative specification of a workflow.

    Attributes:
        name: Human-readable workflow name
        engine: Workflow engine (nextflow | snakemake | wdl)
        entrypoint: Path to workflow script (e.g., main.nf, Snakefile)
        params: Engine-specific parameters
        work_dir: Optional working directory
    """
    name: str
    engine: str
    entrypoint: str
    params: Dict
    work_dir: Optional[str] = None
