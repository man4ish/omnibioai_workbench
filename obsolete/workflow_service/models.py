"""
job/models.py

Defines the Job model for tracking computational or bioinformatics workflow executions.

A Job represents a single workflow run submitted by a user, including its status,
progress, result location, and execution logs. This model is used to monitor and
manage pipelines in the application.

Attributes:
    STATUS_CHOICES (list of tuple): Possible states of a job - queued, running, completed, failed.
"""

from django.db import models
from django.contrib.auth.models import User

class Job(models.Model):
    STATUS_CHOICES = [
        ('queued', 'Queued'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed')
    ]

    plugin_name = models.CharField(max_length=100)
    job_id = models.CharField(max_length=100, unique=True)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='queued')
    progress = models.PositiveIntegerField(default=0)  # 0-100%
    result_path = models.TextField(blank=True, null=True)
    logs = models.TextField(blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return f"{self.plugin_name} - {self.job_id} ({self.status})"
