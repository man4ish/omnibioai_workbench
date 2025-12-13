from django.shortcuts import render
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from omnibioai.services.workflow_service.dashboard.models import WorkflowRun
from .serializers import WorkflowRunSerializer

# View to serve the main dashboard page (React app or static HTML)
def dashboard_view(request):
    return render(request, 'dashboard/index.html')  # React app or static page

# API view to list all workflows
class WorkflowListView(APIView):
    def get(self, request):
        workflows = WorkflowRun.objects.all()
        serializer = WorkflowRunSerializer(workflows, many=True)
        return Response(serializer.data, status=status.HTTP_200_OK)

# API view to get details of a specific workflow
class WorkflowDetailView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            serializer = WorkflowRunSerializer(workflow)
            return Response(serializer.data, status=status.HTTP_200_OK)
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Not found."}, status=status.HTTP_404_NOT_FOUND)

# API view to get progress of a workflow
class WorkflowProgressView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            return Response({"progress": workflow.progress}, status=status.HTTP_200_OK)
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Not found."}, status=status.HTTP_404_NOT_FOUND)

# API view to get logs of a workflow
class WorkflowLogsView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            logs = workflow.logs.splitlines()  # Assuming logs are stored as a long text field
            return Response({"logs": logs}, status=status.HTTP_200_OK)
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Not found."}, status=status.HTTP_404_NOT_FOUND)

# API view to get outputs of a workflow
class WorkflowOutputsView(APIView):
    def get(self, request, workflow_id):
        try:
            workflow = WorkflowRun.objects.get(id=workflow_id)
            outputs = workflow.outputs  # Assuming outputs are stored as JSON field
            return Response({"outputs": outputs}, status=status.HTTP_200_OK)
        except WorkflowRun.DoesNotExist:
            return Response({"detail": "Not found."}, status=status.HTTP_404_NOT_FOUND)
