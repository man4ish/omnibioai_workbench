
from django.urls import path
from . import views

app_name = "gene_annotation"

urlpatterns = [
    path("", views.annotate_view, name="annotate"),
]
