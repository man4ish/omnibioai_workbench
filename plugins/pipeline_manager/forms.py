from django import forms
from .models import Job

class JobForm(forms.ModelForm):
    class Meta:
        model = Job
        fields = ['pipeline_type', 'reference_genome', 'input_file_r1', 'input_file_r2']

    def clean(self):
        cleaned_data = super().clean()
        r1 = cleaned_data.get("input_file_r1")
        r2 = cleaned_data.get("input_file_r2")

        # If pipeline expects paired-end data, ensure both files are provided
        if cleaned_data.get('pipeline_type') == 'rnaseq' and not r2:
            raise forms.ValidationError("Paired-end RNA-seq requires both R1 and R2 FASTQ files.")
