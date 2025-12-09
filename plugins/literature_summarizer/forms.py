from django import forms

class LiteratureForm(forms.Form):
    query = forms.CharField(label="Gene / Variant / Keyword", max_length=100)
    num_results = forms.IntegerField(initial=5, label="Number of Articles")
