from django.urls import path
from . import views
app_name="network_analysis"
urlpatterns = [
    path("", views.upload_network, name="network_analysis_home"),  # handles /network_analysis/
    path("upload/", views.upload_network, name="upload_network"),
]