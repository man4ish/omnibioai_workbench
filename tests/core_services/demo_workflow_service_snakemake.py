"""
demo_workflow_service_snakemake.py

Demonstration script for the OmnibioAI workflow service with Snakemake.

This demo shows how to:
- Define a workflow specification
- Launch a workflow using PipelineManager
- Stream execution state, progress, and logs
- Run without Django, Celery, or external services

Prerequisites:
- Snakemake installed and available in PATH
- Python >= 3.10
"""

from pathlib import Path
from omnibioai.services.workflow_service.pipeline_manager import PipelineManager
from omnibioai.services.workflow_service.models import WorkflowSpec
from omnibioai.services.workflow_service.state import WorkflowState


def ensure_example_workflow() -> Path:
    """
    Create a minimal Snakemake workflow for demo purposes.

    Returns
    -------
    Path
        Absolute path to the Snakefile.
    """
    examples_dir = Path("examples/hello_snakemake")
    examples_dir.mkdir(parents=True, exist_ok=True)

    snakefile_path = examples_dir / "Snakefile"

    if not snakefile_path.exists():
        snakefile_path.write_text(
            """rule all:
    output:
        "hello.txt"
    shell:
        "echo 'Hello from OmnibioAI Snakemake workflow' > hello.txt"
"""
        )

    return snakefile_path.resolve()


def run_demo(mock: bool = False) -> None:
    """
    Run the Snakemake workflow service demo.

    Parameters
    ----------
    mock : bool
        If True, use a mock adapter that doesn't require Snakemake.
    """
    snakefile = ensure_example_workflow()

    spec = WorkflowSpec(
        name="demo_snakemake_workflow",
        engine="snakemake",
        entrypoint=str(snakefile),
        params={},
        work_dir=str(Path("work/demo_snakemake").resolve()),
    )

    manager = PipelineManager()

    # Optional: mock adapter for CI/testing
    if mock:
        from omnibioai.services.workflow_service.adapters.base import BaseWorkflowAdapter

        class MockSnakemakeAdapter(BaseWorkflowAdapter):
            def run(self, entrypoint, work_dir, params):
                for i in range(0, 101, 20):
                    yield i, f"Mock step {i}%"

        manager.executor.ADAPTERS["snakemake"] = lambda: MockSnakemakeAdapter()

    print("\n=== OmnibioAI Workflow Service Demo (Snakemake) ===\n")

    for state, progress, log in manager.launch(spec):
        state_label = state.value if isinstance(state, WorkflowState) else str(state)
        print(f"[{state_label:>9}] {progress:3d}% | {log}")

    print("\n=== Demo completed ===\n")


if __name__ == "__main__":
    # Set mock=True if Snakemake is not installed
    run_demo(mock=False)
