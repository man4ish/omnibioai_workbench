from django.urls import path
from . import views

app_name = 'pipeline_manager'

urlpatterns = [
    path('', views.pipeline_home, name='pipeline_home'),
    path('job/<int:job_id>/', views.job_detail, name='job_detail'),
]
