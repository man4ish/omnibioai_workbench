import pandas as pd
from django.shortcuts import render
from .forms import AnalysisUploadForm

def analysis_upload_view(request):
    context = {}
    if request.method == 'POST':
        form = AnalysisUploadForm(request.POST, request.FILES)
        if form.is_valid():
            analysis_type = form.cleaned_data['analysis_type']
            data_file = request.FILES['data_file']
            
            # Read CSV file to DataFrame
            df = pd.read_csv(data_file)
            
            if analysis_type == 'rna_seq':
                # RNA-Seq dummy logic: classify by mean expression
                mean_expr = df.mean(axis=1, numeric_only=True)
                predictions = ['Tumor' if x > mean_expr.median() else 'Normal' for x in mean_expr]
                context['predictions'] = zip(df.index, predictions)
                context['message'] = "Demo RNA-Seq classification done!"
                
            elif analysis_type == 'dna':
                # DNA dummy logic: classify by median of sum of numeric values per row
                row_sums = df.select_dtypes(include='number').sum(axis=1)
                predictions = ['Mutated' if x > row_sums.median() else 'WildType' for x in row_sums]
                context['predictions'] = zip(df.index, predictions)
                context['message'] = "Demo DNA mutation classification done!"
                
            elif analysis_type == 'geneexpression':
                # Gene expression dummy logic: classify by variance of expression per row
                variances = df.var(axis=1, numeric_only=True)
                predictions = ['HighVar' if x > variances.median() else 'LowVar' for x in variances]
                context['predictions'] = zip(df.index, predictions)
                context['message'] = "Demo Gene Expression variance classification done!"
                
            elif analysis_type == 'dna+rna':
                # Combined DNA+RNA dummy logic: sum of normalized sums (just a dummy approach)
                dna_cols = [col for col in df.columns if 'dna' in col.lower()]
                rna_cols = [col for col in df.columns if 'rna' in col.lower()]
                
                # Sum DNA and RNA numeric values separately
                dna_sum = df[dna_cols].sum(axis=1) if dna_cols else pd.Series(0, index=df.index)
                rna_sum = df[rna_cols].sum(axis=1) if rna_cols else pd.Series(0, index=df.index)
                
                combined_score = dna_sum + rna_sum
                predictions = ['HighScore' if x > combined_score.median() else 'LowScore' for x in combined_score]
                context['predictions'] = zip(df.index, predictions)
                context['message'] = "Demo DNA+RNA combined classification done!"
                
            else:
                context['message'] = f"Received {analysis_type} data file. Prediction logic not implemented yet."
            
            return render(request, 'ml_predictor/results.html', context)
    else:
        form = AnalysisUploadForm()
    return render(request, 'ml_predictor/upload.html', {'form': form})


