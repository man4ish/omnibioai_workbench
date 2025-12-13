"""
state.py

Defines lifecycle states for workflows in the OmnibioAI workflow service.

This module provides the `WorkflowState` enumeration, which represents
the different stages of a workflow's lifecycle from creation to completion.
It is used throughout the workflow management system to track and report
the current state of a workflow.

Key Components:
- `WorkflowState`: Enum class capturing the primary states a workflow can
  be in during its execution.

States:
- `CREATED`: Workflow has been defined but not yet submitted for execution.
- `QUEUED`: Workflow is waiting in the queue to be executed.
- `RUNNING`: Workflow is currently executing.
- `COMPLETED`: Workflow has finished successfully.
- `FAILED`: Workflow execution terminated due to an error.

Design Goals:
- Standardized state management: Provides a consistent way to represent
  workflow states across components (PipelineManager, WorkflowExecutor,
  WorkflowRunner, adapters).
- Readability: Enum values are human-readable and suitable for logging,
  progress reporting, and monitoring.
- Extensibility: New states can be added easily if additional lifecycle
  stages are required.

Example Usage:
    from omnibioai.services.workflow_service.state import WorkflowState

    state = WorkflowState.CREATED
    if state == WorkflowState.CREATED:
        print("Workflow is ready to be queued for execution.")
"""


from enum import Enum


class WorkflowState(Enum):
    CREATED = "created"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
