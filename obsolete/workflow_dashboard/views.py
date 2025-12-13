from django.shortcuts import render, redirect
from plugins.workflow_dashboard.models import WorkflowJob
from omnibioai.services.workflow_service import PipelineManager

pipeline_manager = PipelineManager()

def dashboard(request):
    jobs = WorkflowJob.objects.all().order_by('-created_at')
    return render(request, "workflow_dashboard/dashboard.html", {"jobs": jobs})

def launch_workflow(request):
    if request.method == "POST":
        workflow_type = request.POST.get("workflow_type")
        script_path = request.FILES.get("script_file")  # Save uploaded file
        # Save file locally
        with open(f"/tmp/{script_path.name}", "wb") as f:
            for chunk in script_path.chunks():
                f.write(chunk)

        # Launch workflow
        if workflow_type == "nextflow":
            task_id = pipeline_manager.launch_nextflow(f"/tmp/{script_path.name}")
        elif workflow_type == "snakemake":
            task_id = pipeline_manager.launch_snakemake(f"/tmp/{script_path.name}")
        elif workflow_type == "wdl":
            task_id = pipeline_manager.launch_wdl(f"/tmp/{script_path.name}")

        # Track job
        WorkflowJob.objects.create(task_id=task_id, workflow_type=workflow_type)
        return redirect("workflow_dashboard:dashboard")
