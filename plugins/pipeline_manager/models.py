from django.db import models

class Job(models.Model):
    PIPELINE_CHOICES = [
        ('rnaseq', 'RNA-seq'),
        ('wgs', 'Whole Genome Sequencing'),
        ('wes', 'Whole Exome Sequencing'),
        ('methylation', 'Methylation'),
    ]

    pipeline_type = models.CharField(max_length=20, choices=PIPELINE_CHOICES)
    reference_genome = models.CharField(max_length=100, default='hg38')
    input_file_r1 = models.FileField(upload_to='input_files/')
    input_file_r2 = models.FileField(upload_to='input_files/', blank=True, null=True)  # paired read optional for single-end
    out_dir = models.TextField(blank=True, null=True)
    status = models.CharField(max_length=50, default='PENDING')
    log = models.TextField(blank=True, default='')
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"Job {self.id} ({self.pipeline_type}) - {self.status}"

