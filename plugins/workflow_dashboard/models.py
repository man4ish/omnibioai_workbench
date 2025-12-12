from django.db import models

class WorkflowJob(models.Model):
    task_id = models.CharField(max_length=50, unique=True)
    workflow_type = models.CharField(max_length=20)  # nextflow, snakemake, wdl
    status = models.CharField(max_length=20, default="queued")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    logs = models.TextField(blank=True, null=True)

    def __str__(self):
        return f"{self.workflow_type} - {self.task_id}"
