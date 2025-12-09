from django.urls import path
from . import views

app_name = 'pathway_enrichment'

urlpatterns = [
    path('', views.home, name='upload'),  # 👈 Use 'upload' as name
    path('enrichment/', views.run_enrichr, name='run_enrichment'),
]
