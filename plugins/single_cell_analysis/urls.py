from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from .views import gene_umap, gene_expression  # Import all needed views
from . import views  # Import the views module

app_name = "single_cell_analysis"

urlpatterns = [
    path('', views.upload_file, name='upload_h5ad'),  # Default route → upload
    path('umap/', views.umap_plot, name='umap_plot'),
    path('download/', views.download_metadata, name='download_metadata'),
    path('gene_umap/', gene_umap, name='gene_umap'),
    path('gene_expression/', gene_expression, name='gene_expression'),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)