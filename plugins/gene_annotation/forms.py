from django import forms

class GeneAnnotationForm(forms.Form):
    gene_file = forms.FileField(label="Upload CSV gene expression file")
    prompt = forms.CharField(widget=forms.Textarea, label="Enter prompt for model")
    chart_type = forms.ChoiceField(choices=[('box', 'Box Plot'), ('scatter', 'Scatter Plot'), ('violin', 'Violin Plot'), ('heatmap', 'Heatmap')], initial='box')
    model = forms.ChoiceField(choices=[('llama3', 'LLaMA 3'), ('deepseek', 'DeepSeek')], initial='llama3')
    temperature = forms.FloatField(initial=0.7, min_value=0, max_value=1)
