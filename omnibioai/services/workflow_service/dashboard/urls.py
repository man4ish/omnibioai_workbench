from django.urls import path
from . import views

urlpatterns = [
    # Render the main dashboard view (serving React app or a static dashboard page)
    path('', views.dashboard_view, name='dashboard'),

    # API routes for workflow progress, logs, and outputs
    path('api/workflows/', views.WorkflowListView.as_view(), name='workflow-list'),
    path('api/workflows/<int:workflow_id>/', views.WorkflowDetailView.as_view(), name='workflow-detail'),
    path('api/workflows/<int:workflow_id>/progress/', views.WorkflowProgressView.as_view(), name='workflow-progress'),
    path('api/workflows/<int:workflow_id>/logs/', views.WorkflowLogsView.as_view(), name='workflow-logs'),
    path('api/workflows/<int:workflow_id>/outputs/', views.WorkflowOutputsView.as_view(), name='workflow-outputs'),
]
