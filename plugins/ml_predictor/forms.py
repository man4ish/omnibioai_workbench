from django import forms

ANALYSIS_CHOICES = [
    ('rna_seq', 'RNA-Seq'),
    ('dna', 'DNA'),
    ('geneexpression', 'Gene Expression'),
    ('dna+rna', 'DNA + RNA'),
]

class AnalysisUploadForm(forms.Form):
    analysis_type = forms.ChoiceField(choices=ANALYSIS_CHOICES, label="Select Analysis Type")
    data_file = forms.FileField(label="Upload Data File")
