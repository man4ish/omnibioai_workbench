from django.urls import path
from .views import analysis_upload_view

app_name = 'ml_predictor'

urlpatterns = [
    path('', analysis_upload_view, name='ml_predictor_upload'),  # <-- name here must match
]