"""
URL configuration for omnibioai project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.contrib import admin
from django.urls import path, include
from django.views.generic import RedirectView

urlpatterns = [
    path("admin/", admin.site.urls),
    path('plugins/igv_viewer/', include('plugins.igv_viewer.urls')),  # <-- add this
    path('plugins/pipeline_manager/', include('plugins.pipeline_manager.urls')),  # add this
    path('plugins/gene_annotation/', include('plugins.gene_annotation.urls')),
    path('plugins/ml_predictor/', include('plugins.ml_predictor.urls')),
    path('plugins/variant_annotation/', include('plugins.variant_annotation.urls')),
    path('plugins/pathway_enrichment/', include('plugins.pathway_enrichment.urls')),
    path('plugins/literature_summarizer/', include('plugins.literature_summarizer.urls')),
    path('plugins/network_analysis/', include('plugins.network_analysis.urls')),
    path('plugins/single_cell_analysis/', include('plugins.single_cell_analysis.urls')),

    path('', include('plugins.home.urls')),   # root path
]
