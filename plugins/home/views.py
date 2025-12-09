from django.shortcuts import render

def home_view(request):
    apps = [
        {'name': 'ML Predictor', 'url': '/ml_predictor/'},
        {'name': 'Pipeline Manager', 'url': '/pipeline_manager/'},
        {'name': 'IGV Viewer', 'url': '/igv_viewer/'},
        # Add more apps here as you create them
    ]
    return render(request, 'home/home.html', {'apps': apps})