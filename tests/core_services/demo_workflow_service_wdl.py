"""
demo_workflow_service_wdl.py

Demonstration script for the OmnibioAI workflow service with WDL (Cromwell).

This demo shows how to:
- Define a workflow specification
- Launch a workflow using PipelineManager
- Stream execution state, progress, and logs
- Run without Django, Celery, or external services

Prerequisites:
- Java available in PATH
- Python >= 3.10
"""

from pathlib import Path
from omnibioai.services.workflow_service.pipeline_manager import PipelineManager
from omnibioai.services.workflow_service.models import WorkflowSpec
from omnibioai.services.workflow_service.state import WorkflowState


def ensure_example_workflow() -> Path:
    """
    Create a minimal WDL workflow for demo purposes.

    Returns
    -------
    Path
        Absolute path to the WDL workflow file.
    """
    examples_dir = Path("examples/hello_wdl")
    examples_dir.mkdir(parents=True, exist_ok=True)

    wdl_file = examples_dir / "hello.wdl"

    if not wdl_file.exists():
        wdl_file.write_text(
            """
version 1.0

workflow hello_wdl {
    call hello
}

task hello {
    command {
        echo "Hello from OmnibioAI WDL workflow"
    }
    output {
        String message = read_string(stdout())
    }
}
"""
        )

    return wdl_file.resolve()


def run_demo() -> None:
    """
    Run the WDL workflow service demo.
    Automatically falls back to a mock adapter if Cromwell is missing.
    """
    wdl_file = ensure_example_workflow()
    cromwell_jar = Path("bin/cromwell.jar").resolve()
    use_mock = not cromwell_jar.exists()

    spec = WorkflowSpec(
        name="demo_wdl_workflow",
        engine="wdl",
        entrypoint=str(wdl_file),
        params={},
        work_dir=str(Path("work/demo_wdl").resolve()),
    )

    manager = PipelineManager()

    if use_mock:
        from omnibioai.services.workflow_service.adapters.base import BaseWorkflowAdapter

        class MockWDLAdapter(BaseWorkflowAdapter):
            def run(self, entrypoint, work_dir, params):
                for i in range(0, 101, 20):
                    yield i, f"Mock step {i}%"

        manager.executor.ADAPTERS["wdl"] = lambda: MockWDLAdapter()
        print("Cromwell JAR not found; using mock WDL adapter.\n")

    else:
        from omnibioai.services.workflow_service.adapters.wdl_adapter import WDLAdapter
        manager.executor.ADAPTERS["wdl"] = lambda: WDLAdapter(
            cromwell_jar=str(cromwell_jar),  # Correct argument name
            java_bin="java"
        )

    print("\n=== OmnibioAI Workflow Service Demo (WDL) ===\n")

    for state, progress, log in manager.launch(spec):
        state_label = state.value if isinstance(state, WorkflowState) else str(state)
        print(f"[{state_label:>9}] {progress:3d}% | {log}")

    print("\n=== Demo completed ===\n")


if __name__ == "__main__":
    run_demo()
