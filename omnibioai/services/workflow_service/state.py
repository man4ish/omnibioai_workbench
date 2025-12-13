"""
state.py

Defines lifecycle states for workflow execution.
"""

from enum import Enum


class WorkflowState(Enum):
    CREATED = "created"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
