from rest_framework import serializers
from omnibioai.services.workflow_service.models import WorkflowRun

class WorkflowRunSerializer(serializers.ModelSerializer):
    """
    Serializer for the WorkflowRun model.
    """
    class Meta:
        model = WorkflowRun
        fields = ['id', 'workflow_name', 'engine', 'entrypoint', 'parameters', 
                  'start_time', 'end_time', 'state', 'outputs', 'logs']
