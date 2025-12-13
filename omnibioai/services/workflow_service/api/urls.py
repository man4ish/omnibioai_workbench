from django.urls import path
from . import views

urlpatterns = [
    # List all workflows (GET)
    path('workflows/', views.WorkflowListView.as_view(), name='workflow-list'),
    
    # Get details for a specific workflow (GET)
    path('workflows/<int:workflow_id>/', views.WorkflowDetailView.as_view(), name='workflow-detail'),
    
    # Get logs for a specific workflow (GET)
    path('workflows/<int:workflow_id>/logs/', views.WorkflowLogsView.as_view(), name='workflow-logs'),
    
    # Get outputs for a specific workflow (GET)
    path('workflows/<int:workflow_id>/outputs/', views.WorkflowOutputsView.as_view(), name='workflow-outputs'),
    
    # Get progress for a specific workflow (GET)
    path('workflows/<int:workflow_id>/progress/', views.WorkflowProgressView.as_view(), name='workflow-progress'),
]
