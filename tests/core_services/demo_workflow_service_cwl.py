"""
Demo: OmnibioAI Workflow Service (CWL)

This demo executes a CWL workflow using PipelineManager.
It automatically creates a minimal CWL file if missing and
uses an absolute path to avoid file-not-found issues.
"""

import os
from omnibioai.services.workflow_service.pipeline_manager import PipelineManager
from omnibioai.services.workflow_service.models import WorkflowSpec
from omnibioai.services.workflow_service.state import WorkflowState

# ------------------------------
# Ensure CWL workflow exists
# ------------------------------
workflow_dir = "tests/data"
workflow_file = "hello_world.cwl"
workflow_path = os.path.join(workflow_dir, workflow_file)
workflow_path = os.path.abspath(workflow_path)

os.makedirs(workflow_dir, exist_ok=True)

if not os.path.exists(workflow_path):
    print(f"[INFO] Creating minimal CWL workflow at {workflow_path}")
    with open(workflow_path, "w") as f:
        f.write(
            """cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo

inputs:
  message:
    type: string
    inputBinding:
      position: 1

outputs: {}
"""
        )

# ------------------------------
# Set up PipelineManager and WorkflowSpec
# ------------------------------
manager = PipelineManager()

spec = WorkflowSpec(
    name="CWL_HelloWorld",
    engine="cwl",
    entrypoint=workflow_path,  # absolute path
    params={"message": "Hello from OmnibioAI CWL adapter"},
    work_dir="./work/cwl_demo"
)

# ------------------------------
# Run workflow
# ------------------------------
print("=== OmnibioAI Workflow Service Demo (CWL) ===")

try:
    for state, progress, log in manager.launch(spec):
        print(f"[{state.name:10s}] {progress:3d}% | {log}")

        if state in (WorkflowState.FAILED, WorkflowState.COMPLETED):
            break

    print("=== CWL Workflow Finished ===")

except Exception as e:
    print("=== CWL Workflow Failed ===")
    print(str(e))
