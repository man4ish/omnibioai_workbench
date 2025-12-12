# workflow_service/__init__.py
from .pipeline_manager import PipelineManager
from .executor import app as celery_app

__all__ = ["PipelineManager", "celery_app"]
