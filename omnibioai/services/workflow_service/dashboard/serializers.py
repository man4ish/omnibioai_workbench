from rest_framework import serializers
from omnibioai.services.workflow_service.dashboard.models import WorkflowRun

class WorkflowRunSerializer(serializers.ModelSerializer):
    class Meta:
        model = WorkflowRun
        fields = '__all__'
