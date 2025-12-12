from django.urls import path
from . import views

app_name = "workflow_dashboard"

urlpatterns = [
    path("", views.dashboard, name="dashboard"),
    path("launch/", views.launch_workflow, name="launch_workflow"),
]
