"""
models.py

Defines core data models for workflow execution.
These models are framework-agnostic and contain no execution logic.
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
