from django import forms

class UploadH5ADForm(forms.Form):
    file = forms.FileField(label='Upload .h5ad file')

class GeneExpressionForm(forms.Form):
    gene = forms.CharField(label="Gene name", max_length=100)


class SelectColorForm(forms.Form):
    color_by = forms.ChoiceField(
        label="Color by",
        choices=[],
        required=False,
        help_text="Select a column to color the plot by"
    )
    method = forms.ChoiceField(
        label="Dimensionality Reduction Method",
        choices=[('umap', 'UMAP'), ('tsne', 't-SNE'), ('pca', 'PCA')],
        required=True,
        initial='umap'
    )

    def __init__(self, *args, **kwargs):
        color_choices = kwargs.pop('color_choices', [])
        super().__init__(*args, **kwargs)

        # Add a default empty choice to color_by to allow no selection
        choices = [('', 'Select color (optional)')] + [(c, c) for c in color_choices]
        self.fields['color_by'].choices = choices