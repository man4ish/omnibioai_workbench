from django.http import JsonResponse
from django.shortcuts import render

# Genome list data - you can extend this later or load from DB/file
GENOME_LIST = [
    {"id": "hg19", "name": "Human (hg19)"},
    {"id": "hg38", "name": "Human (hg38)"},
    {"id": "mm10", "name": "Mouse (mm10)"},
]

def igv_home(request):
    # Just render the page; genome list is fetched by AJAX
    return render(request, "igv_viewer/igv_home.html")

def genome_list(request):
    return JsonResponse(GENOME_LIST, safe=False)
