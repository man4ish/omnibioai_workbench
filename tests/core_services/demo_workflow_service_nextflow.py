"""
demo_workflow_service_nextflow.py

Demonstration script for the OmnibioAI workflow service.

This demo shows how to:
- Define a workflow specification
- Launch a workflow using PipelineManager
- Stream execution state, progress, and logs
- Run without Django, Celery, or external services

Prerequisites:
- Nextflow installed and available in PATH
- Python >= 3.10 (tested on 3.13)
"""

from pathlib import Path
from omnibioai.services.workflow_service.pipeline_manager import PipelineManager
from omnibioai.services.workflow_service.models import WorkflowSpec
from omnibioai.services.workflow_service.state import WorkflowState


def ensure_example_workflow() -> Path:
    """
    Create a minimal Nextflow workflow for demo purposes.

    Returns
    -------
    Path
        Absolute path to the workflow file.
    """
    examples_dir = Path("examples")
    examples_dir.mkdir(parents=True, exist_ok=True)

    workflow_path = examples_dir / "hello_nextflow.nf"

    if not workflow_path.exists():
        workflow_path.write_text(
            """
process hello {
    output:
    stdout

    script:
    '''
    echo "Hello from OmnibioAI workflow service"
    '''
}

workflow {
    hello()
}
"""
        )

    return workflow_path.resolve()  # return absolute path


def run_demo(mock: bool = False) -> None:
    """
    Run the workflow service demo.

    Parameters
    ----------
    mock : bool
        If True, use a mock adapter that doesn't require Nextflow.
    """
    workflow_file = ensure_example_workflow()

    spec = WorkflowSpec(
        name="demo_nextflow_workflow",
        engine="nextflow",
        entrypoint=str(workflow_file),
        params={},
        work_dir=str(Path("work/demo_nextflow").resolve()),
    )

    manager = PipelineManager()

    # Optional: swap in mock adapter
    if mock:
        from omnibioai.services.workflow_service.adapters.base import BaseWorkflowAdapter

        class MockNextflowAdapter(BaseWorkflowAdapter):
            def run(self, entrypoint, work_dir, params):
                for i in range(0, 101, 20):
                    yield i, f"Mock step {i}%"

        manager.executor.ADAPTERS["nextflow"] = lambda: MockNextflowAdapter()

    print("\n=== OmnibioAI Workflow Service Demo (Nextflow) ===\n")

    for state, progress, log in manager.launch(spec):
        state_label = (
            state.value if isinstance(state, WorkflowState) else str(state)
        )
        print(f"[{state_label:>9}] {progress:3d}% | {log}")

    print("\n=== Demo completed ===\n")


if __name__ == "__main__":
    # Pass mock=True if Nextflow is not installed
    run_demo(mock=False)
