"""
workflow_executor.py

Defines Celery tasks and utilities for executing bioinformatics workflows
(Nextflow, Snakemake, and WDL) asynchronously, tracking progress in the
database, and sending real-time updates via WebSockets.

This module provides:
- Celery app configuration for distributed task execution.
- Workflow adapters for Nextflow, Snakemake, and WDL pipelines.
- Utility functions to create Job entries in the database.
- Functions to update Job progress and notify connected clients via Channels.
- Celery tasks to run workflows asynchronously and update Job status/logs.

Key Components:
- `app` (Celery): Configured Celery application with Redis backend.
- `create_job(plugin_name, user_id) -> Job`: Create a Job entry in the database.
- `notify_progress(job)`: Send WebSocket updates about the Job.
- `update_job_progress(job, progress, logs)`: Update Job progress and logs.
- `run_nextflow_workflow`, `run_snakemake_workflow`, `run_wdl_workflow`:
  Celery tasks that execute workflows via their respective adapters.

Jobs are tracked using the `Job` Django model, with status transitions
from 'queued' → 'running' → 'completed'/'failed'.
"""

import uuid
from celery import Celery, current_task
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from django.contrib.auth.models import User

from .adapters import NextflowAdapter, SnakemakeAdapter, WDLAdapter
from .models import Job

# -----------------------------
# Celery app
# -----------------------------
app = Celery(
    "workflow_executor",
    broker="redis://localhost:6379/0",
    backend="redis://localhost:6379/0"
)

# -----------------------------
# Initialize workflow adapters
# -----------------------------
nextflow_adapter = NextflowAdapter()
snakemake_adapter = SnakemakeAdapter()
wdl_adapter = WDLAdapter()

# -----------------------------
# Utility: Create Job in DB
# -----------------------------
def create_job(plugin_name: str, user_id: int) -> Job:
    user = User.objects.get(id=user_id)
    job = Job.objects.create(
        plugin_name=plugin_name,
        job_id=str(uuid.uuid4()),
        user=user,
        status="queued",
        progress=0
    )
    return job

# -----------------------------
# Utility: Send WebSocket updates
# -----------------------------
def notify_progress(job: Job):
    channel_layer = get_channel_layer()
    async_to_sync(channel_layer.group_send)(
        "workflow_updates",
        {
            "type": "workflow_update",
            "data": {
                "task_id": job.job_id,
                "plugin_name": job.plugin_name,
                "status": job.status,
                "progress": job.progress,
                "logs": job.logs or "",
            },
        },
    )

# -----------------------------
# Utility: Update job progress
# -----------------------------
def update_job_progress(job: Job, progress: int, logs: str = None):
    job.progress = progress
    if logs:
        job.logs = (job.logs or "") + f"\n{logs}"
    if progress >= 100:
        job.status = "completed"
    elif progress > 0:
        job.status = "running"
    job.save()
    notify_progress(job)

# -----------------------------
# Celery tasks for workflows
# -----------------------------
@app.task(bind=True)
def run_nextflow_workflow(self, script_path, user_id, work_dir=None, params=None):
    job = create_job("Nextflow", user_id)
    try:
        for step, log in nextflow_adapter.run_workflow(script_path, work_dir, params, progress=True):
            update_job_progress(job, step, logs=log)
    except Exception as e:
        job.status = "failed"
        job.logs = (job.logs or "") + f"\n{str(e)}"
        job.save()
        raise

@app.task(bind=True)
def run_snakemake_workflow(self, snakefile_path, user_id, work_dir=None, config_file=None):
    job = create_job("Snakemake", user_id)
    try:
        for step, log in snakemake_adapter.run_workflow(snakefile_path, work_dir, config_file, progress=True):
            update_job_progress(job, step, logs=log)
    except Exception as e:
        job.status = "failed"
        job.logs = (job.logs or "") + f"\n{str(e)}"
        job.save()
        raise

@app.task(bind=True)
def run_wdl_workflow(self, wdl_file, user_id, inputs_json=None, options_json=None, work_dir=None):
    job = create_job("WDL", user_id)
    try:
        for step, log in wdl_adapter.run_workflow(wdl_file, inputs_json, options_json, work_dir, progress=True):
            update_job_progress(job, step, logs=log)
    except Exception as e:
        job.status = "failed"
        job.logs = (job.logs or "") + f"\n{str(e)}"
        job.save()
        raise
