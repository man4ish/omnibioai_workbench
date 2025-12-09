from django import forms

class GeneListUploadForm(forms.Form):
    gene_list = forms.FileField(label="Upload Gene List File (TSV with gene and pvalue)")
