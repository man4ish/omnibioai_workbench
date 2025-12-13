from rest_framework.views import APIView
from rest_framework.response import Response
from omnibioai.services.workflow_service.dashboard.models import WorkflowRun  # Fix the typo here
from omnibioai.services.workflow_service.api.serializers import WorkflowRunSerializer

# List all workflows
class WorkflowListView(APIView):
    def get(self, request):
        workflows = WorkflowRun.objects.all()
        serializer = WorkflowRunSerializer(workflows, many=True)
        return Response(serializer.data)

# Get details for a specific workflow
class WorkflowDetailView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            serializer = WorkflowRunSerializer(workflow)
            return Response(serializer.data)
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Workflow not found"}, status=404)

# Get logs for a specific workflow
class WorkflowLogsView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            logs = workflow.logs.splitlines()  # assuming logs are stored as a long text field
            return Response({"logs": logs})
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Workflow not found"}, status=404)

# Get outputs for a specific workflow
class WorkflowOutputsView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            outputs = workflow.outputs  # Assuming it's a JSON field with outputs
            return Response({"outputs": outputs})
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Workflow not found"}, status=404)

# Get progress for a specific workflow
class WorkflowProgressView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            return Response({"progress": workflow.progress})
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Workflow not found"}, status=404)
