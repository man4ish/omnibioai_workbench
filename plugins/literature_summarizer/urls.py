from django.urls import path
from . import views

app_name = 'literature_summarizer'  # <-- THIS IS ESSENTIAL

urlpatterns = [
    path("", views.summarize_view, name="literature_summary"), 
    path("trends/", views.trend_analysis, name="trend_analysis"),
    path("cooccurrence/", views.keyword_network, name="keyword_network"), 
    path("annotate/", views.pubtator_annotate, name="pubtator_annotate"),
]
