import os

import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for macOS and server compatibility
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

from django.conf import settings
from django.contrib import messages
from django.http import HttpResponse
from django.shortcuts import render, redirect
from django.utils.safestring import mark_safe
from django.views.decorators.csrf import csrf_exempt

from .forms import UploadH5ADForm, SelectColorForm, GeneExpressionForm



# Define upload directory for h5ad files
UPLOAD_DIR = os.path.join(settings.BASE_DIR, 'data', 'single_cell')
os.makedirs(UPLOAD_DIR, exist_ok=True)

# Define directory for UMAP plot
FIG_DIR = os.path.join(settings.BASE_DIR, 'static', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

import os
import scanpy as sc
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import plotly.express as px

from django.conf import settings
from django.contrib import messages
from django.http import HttpResponse
from django.shortcuts import render, redirect
from django.utils.safestring import mark_safe

from .forms import UploadH5ADForm, SelectColorForm

# Define upload directory
UPLOAD_DIR = os.path.join(settings.BASE_DIR, 'data', 'single_cell')
os.makedirs(UPLOAD_DIR, exist_ok=True)

def upload_file(request):
    if request.method == 'POST':
        form = UploadH5ADForm(request.POST, request.FILES)
        if form.is_valid():
            # saves the uploaded file
            ...
            return redirect('single_cell_analysis:umap_plot')
    else:
        form = UploadH5ADForm()
    return render(request, 'single_cell_analysis/upload.html', {'form': form})



def umap_plot(request):
    file_path = request.session.get("h5ad_file")
    if not file_path or not os.path.exists(file_path):
        messages.error(request, "Please upload an .h5ad first.")
        return redirect("single_cell_analysis:upload_h5ad")

    adata = sc.read_h5ad(file_path)

    # Basic preprocessing
    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata)
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata)

    # Get user selections
    method = request.GET.get("method", "umap")
    colors = list(adata.obs.columns)
    default_color = colors[0] if colors else None
    color_by = request.GET.get("color_by", default_color)

    # Compute chosen embedding if needed
    if method == "umap":
        sc.tl.umap(adata)
        coords = adata.obsm["X_umap"]
        xcol, ycol = "UMAP1", "UMAP2"
    elif method == "tsne":
        sc.tl.tsne(adata)
        coords = adata.obsm["X_tsne"]
        xcol, ycol = "tSNE1", "tSNE2"
    elif method == "pca":
        # PCA coords exist; take only first two PCs
        coords = adata.obsm["X_pca"][:, :2]
        xcol, ycol = "PCA1", "PCA2"
    else:
        return render(request, "single_cell_analysis/umap.html", {
            "error": f"Unknown method '{method}'"
        })

    # Build DataFrame of coordinates
    df = pd.DataFrame(coords, columns=[xcol, ycol], index=adata.obs_names)

    # Only add color column if a valid column name was provided
    if color_by:
        df[color_by] = adata.obs[color_by].values

    # Create Plotly scatter
    fig = px.scatter(
        df,
        x=xcol,
        y=ycol,
        color=color_by if color_by else None,
        title=f"{method.upper()} colored by {color_by or 'none'}",
        hover_name=df.index,
        width=800,
        height=600
    )
    plot_html = fig.to_html(full_html=False)

    # Build form with both method & color fields
    form = SelectColorForm(
        color_choices=colors,
        initial={"method": method, "color_by": color_by}
    )

    return render(request, "single_cell_analysis/umap.html", {
        "form": form,
        "filename": os.path.basename(file_path),
        "plot_html": mark_safe(plot_html),
    })

def download_metadata(request):
    file_path = request.session.get('h5ad_file')
    if not file_path or not os.path.exists(file_path):
        return HttpResponse("No file found.", status=404)

    adata = sc.read_h5ad(file_path)
    df = adata.obs
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename=metadata.csv'
    df.to_csv(path_or_buf=response)
    return response

@csrf_exempt
def gene_umap(request):
    if request.method == "POST":
        gene_name = request.POST.get("gene_name", None)
        if not gene_name:
            return render(request, "single_cell_analysis/gene_umap.html", {"error": "Please enter a gene name."})

        # Your UMAP plotting code here...
        # For demonstration, create a dummy plot
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [4, 5, 6])
        ax.set_title(f"UMAP plot for {gene_name}")

        # Ensure media directory exists
        media_dir = settings.MEDIA_ROOT
        if not os.path.exists(media_dir):
            os.makedirs(media_dir)

        img_filename = "gene_umap.png"
        img_path = os.path.join(media_dir, img_filename)
        fig.savefig(img_path)
        plt.close(fig)  # Close plot to free memory

        # Pass relative URL to template (use MEDIA_URL)
        img_url = os.path.join(settings.MEDIA_URL, img_filename)

        return render(request, "single_cell_analysis/gene_umap.html", {
            "img_path": img_url,
            "gene_name": gene_name,
        })

    # GET request - just render the input form
    return render(request, "single_cell_analysis/gene_umap.html")

def gene_expression(request):
    plot_html = None
    error = None

    if request.method == 'POST':
        h5ad_file = request.FILES.get('h5ad_file')
        gene_name = request.POST.get('gene_name', '').strip()  # <-- use form input here

        if not h5ad_file:
            error = "No file uploaded"
        elif not gene_name:
            error = "Please enter a gene name."
        else:
            try:
                temp_path = os.path.join(settings.MEDIA_ROOT, h5ad_file.name)
                with open(temp_path, 'wb+') as dest:
                    for chunk in h5ad_file.chunks():
                        dest.write(chunk)

                adata = sc.read_h5ad(temp_path)

                if gene_name not in adata.var_names:
                    error = f"Gene {gene_name} not found in dataset."
                else:
                    if 'clusters' not in adata.obs:
                        adata.obs['clusters'] = '0'

                    df = adata.to_df()
                    df['cluster'] = adata.obs['clusters'].values

                    mean_expr = df.groupby('cluster')[gene_name].mean().reset_index()

                    fig = px.bar(mean_expr, x='cluster', y=gene_name,
                                 labels={'cluster': 'Cluster', gene_name: 'Mean Expression'},
                                 title=f'Mean Expression of {gene_name} Across Clusters')

                    plot_html = fig.to_html(full_html=False)

                os.remove(temp_path)

            except Exception as e:
                error = f"Error processing file: {str(e)}"

    return render(request, 'single_cell_analysis/gene_expression.html', {
        'plot_html': plot_html,
        'error': error,
    })

