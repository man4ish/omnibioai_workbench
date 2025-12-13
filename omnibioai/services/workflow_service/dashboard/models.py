from django.db import models

class WorkflowRun(models.Model):
    # Unique ID for each workflow run (useful for tracking)
    run_id = models.CharField(max_length=255, unique=True, null=False, blank=False)

    # Workflow name and engine type
    workflow_name = models.CharField(max_length=255)
    engine = models.CharField(max_length=50)
    
    # Entry point (main file or script of the workflow)
    entrypoint = models.TextField()
    
    # Parameters used in the workflow
    parameters = models.JSONField()

    # Timestamps
    start_time = models.DateTimeField(auto_now_add=True)  # Set automatically when a workflow starts
    end_time = models.DateTimeField(null=True, blank=True)  # Can be set when the workflow ends
    
    # Workflow state (running, completed, failed)
    state = models.CharField(max_length=20)
    
    # Outputs (could be JSON or file paths depending on your data)
    outputs = models.JSONField(null=True, blank=True)
    
    # Logs captured during the workflow execution
    logs = models.TextField()

    # Optional: Store progress for running workflows (0-100%)
    progress = models.IntegerField(default=0, null=False, blank=False)
    
    def __str__(self):
        return f"{self.workflow_name} - {self.state} - {self.run_id}"

    # Optional: If you need to track unique workflow runs
    class Meta:
        indexes = [
            models.Index(fields=['workflow_name']),
            models.Index(fields=['state']),
            models.Index(fields=['start_time']),
        ]
