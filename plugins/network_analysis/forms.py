from django import forms

class UploadNetworkForm(forms.Form):
    file = forms.FileField(label="Upload interaction file (CSV or TSV)")