from django.shortcuts import render, redirect, get_object_or_404
from .models import Job
from .forms import JobForm
from .tasks import run_rnaseq_nextflow_pipeline, run_wgs_pipeline, run_wes_pipeline, run_methylation_pipeline

def pipeline_home(request):
    if request.method == 'POST':
        form = JobForm(request.POST, request.FILES)
        if form.is_valid():
            job = form.save()

            # Trigger Celery task based on pipeline_type
            if job.pipeline_type == 'rnaseq':
                run_rnaseq_nextflow_pipeline.delay(job.id)
            elif job.pipeline_type == 'wgs':
                run_wgs_pipeline.delay(job.id)
            elif job.pipeline_type == 'wes':
                run_wes_pipeline.delay(job.id)
            elif job.pipeline_type == 'methylation':
                run_methylation_pipeline.delay(job.id)

            return redirect('pipeline_manager:pipeline_home')
    else:
        form = JobForm()

    jobs = Job.objects.all().order_by('-created_at')  # or any other ordering

    return render(request, 'pipeline_manager/home.html', {'form': form, 'jobs': jobs})

def job_detail(request, job_id):
    job = get_object_or_404(Job, id=job_id)
    return render(request, 'pipeline_manager/job_detail.html', {'job': job})
