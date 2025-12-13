"""
visualization_service.py

Core plotting and visualization utilities for OmnibioAI.

This module provides a comprehensive suite of plotting functions for genomic,
omics, and general data analysis workflows. It supports both static (Matplotlib)
and interactive (Plotly) outputs, with optional PNG export via Plotly/Kaleido.
Functions handle typical bioinformatics plots (PCA, Volcano, Manhattan, Heatmap,
Boxplot, Lollipop, Bar, Histogram, Pie/Donut, Line/Scatter) and allow reading
from a pandas DataFrame or CSV path.

Features:
---------
- Headless plotting for server environments (Matplotlib 'Agg' backend)
- PCA computation with Plotly and Matplotlib outputs
- Volcano, Manhattan, Heatmap, Boxplot, Bar, Pie/Donut, Line/Scatter, Lollipop plots
- IGV.js-compatible track metadata generation
- Convenience wrapper to select the preferred output format (plotly JSON, PNG bytes)
- Automatic results directory creation and safe file saving
- Optional sample or group annotations for plots
- Return metadata for downstream processing

Key Functions:
--------------
- plot_pca(): Principal Component Analysis plot
- plot_volcano(): Volcano plot for differential analysis
- plot_manhattan(): Manhattan plot for GWAS or genomic data
- plot_heatmap(): Heatmap from numeric matrices
- plot_boxplot(): Grouped or univariate boxplot
- plot_histogram(): Histogram for numeric columns
- plot_bar(): Bar plot for categorical summaries
- plot_pie_donut(): Pie or donut chart
- plot_line_scatter(): Line or scatter plot
- plot_lollipop(): Lollipop plot for category vs value
- igv_track_meta(): Generate IGV.js track metadata
- choose_best_output(): Select preferred output format from plotting result

Usage:
------
>>> from omnibioai.services.visualization_service import plot_pca
>>> import pandas as pd
>>> df = pd.DataFrame(...)  # your data
>>> res = plot_pca(df, sample_column="sample_id", save_path_prefix="pca_demo")
>>> print(res["explained_variance_ratio"])
>>> if "file_path" in res:
>>>     print("Saved plot to:", res["file_path"])
"""


import os
import io
import json
from typing import Optional, Dict, Any, Tuple, Union, List

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use("Agg")  # headless backend for servers
import matplotlib.pyplot as plt

import plotly.graph_objects as go
import plotly.express as px

# Optional: Plotly PNG export requires kaleido
_PLOTLY_PNG_SUPPORTED = True
try:
    # Plotly uses kaleido under the hood; this import is just a check
    import kaleido  # noqa: F401
except Exception:
    _PLOTLY_PNG_SUPPORTED = False

# Default results directory
DEFAULT_RESULTS_DIR = os.environ.get("OMNIBIOAI_RESULTS_DIR", os.path.join("data", "results"))
os.makedirs(DEFAULT_RESULTS_DIR, exist_ok=True)


# -------------------------
# Utilities
# -------------------------
def _read_df(df: Optional[pd.DataFrame], csv_path: Optional[str]) -> pd.DataFrame:
    if df is not None:
        if not isinstance(df, pd.DataFrame):
            raise ValueError("`df` must be a pandas DataFrame")
        return df.copy()
    if csv_path is not None:
        return pd.read_csv(csv_path)
    raise ValueError("Either `df` or `csv_path` must be provided.")


def _save_bytes_to_file(b: bytes, prefix: str, ext: str = "png", results_dir: str = DEFAULT_RESULTS_DIR) -> str:
    os.makedirs(results_dir, exist_ok=True)
    fname = f"{prefix}.{ext}"
    path = os.path.join(results_dir, fname)
    # ensure unique if exists
    base, counter = fname, 1
    while os.path.exists(path):
        fname = f"{prefix}_{counter}.{ext}"
        path = os.path.join(results_dir, fname)
        counter += 1
    with open(path, "wb") as fh:
        fh.write(b)
    return path


def _fig_to_png_bytes_matplotlib(fig: matplotlib.figure.Figure, dpi: int = 150) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def _plotly_fig_to_png_bytes(fig: go.Figure, format: str = "png") -> Optional[bytes]:
    if not _PLOTLY_PNG_SUPPORTED:
        return None
    try:
        # to_image returns bytes
        img = fig.to_image(format=format)
        return img
    except Exception:
        return None


def _fig_to_plotly_json(fig: go.Figure) -> Dict[str, Any]:
    return fig.to_dict()


# -------------------------
# Core plotting functions
# -------------------------
def plot_pca(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    n_components: int = 2,
    center: bool = True,
    scale: bool = True,
    sample_column: Optional[str] = None,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    """
    Compute PCA and return Plotly JSON and PNG.

    Returns dict:
      - plotly: plotly figure dict
      - png_bytes: matplotlib PNG bytes
      - plotly_png_bytes: plotly-exported PNG bytes (None if kaleido missing)
      - explained_variance_ratio: list
      - file_path: saved PNG path if save_path_prefix provided
    """
    return_types = return_types or ["plotly", "png_bytes"]
    df = _read_df(df, csv_path)

    # handle sample column
    sample_names = None
    if sample_column and sample_column in df.columns:
        sample_names = df[sample_column].astype(str).tolist()
        df = df.drop(columns=[sample_column])

    numeric_df = df.select_dtypes(include=[np.number])
    if numeric_df.shape[1] == 0:
        raise ValueError("No numeric columns found for PCA.")

    X = numeric_df.values.astype(float)
    if center:
        X = X - np.nanmean(X, axis=0)
    if scale:
        s = np.nanstd(X, axis=0)
        s[s == 0] = 1.0
        X = X / s

    pca = PCA(n_components=n_components)
    comps = pca.fit_transform(np.nan_to_num(X))
    comp_cols = [f"PC{i+1}" for i in range(comps.shape[1])]
    df_pca = pd.DataFrame(comps, columns=comp_cols)
    if sample_names is not None:
        df_pca['sample'] = sample_names

    # Plotly interactive
    fig_plotly = None
    if "plotly" in return_types:
        if comps.shape[1] >= 2:
            fig_plotly = px.scatter(
                df_pca,
                x='PC1',
                y='PC2',
                hover_name='sample' if 'sample' in df_pca.columns else None,
                title='PCA: PC1 vs PC2'
            )
        else:
            fig_plotly = go.Figure()
            fig_plotly.add_trace(go.Scatter(x=df_pca['PC1'], y=[0]*len(df_pca), mode='markers'))
            fig_plotly.update_layout(title='PCA: PC1')

    # Matplotlib static PNG
    fig_mpl = plt.figure(figsize=(6, 4))
    ax = fig_mpl.add_subplot(111)
    if comps.shape[1] >= 2:
        ax.scatter(comps[:, 0], comps[:, 1], s=20)
        if sample_names:
            for i, txt in enumerate(sample_names):
                ax.annotate(txt, (comps[i, 0], comps[i, 1]), fontsize=6)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("PCA: PC1 vs PC2")
    else:
        ax.scatter(comps[:, 0], np.zeros(len(comps)))
        ax.set_xlabel("PC1")
        ax.set_title("PCA: PC1")

    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)

    # Try Plotly PNG
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly) if fig_plotly is not None else None

    out = {
        "plotly": _fig_to_plotly_json(fig_plotly) if fig_plotly is not None else None,
        "png_bytes": png_bytes,
        "plotly_png_bytes": plotly_png,
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "meta": {"n_components": n_components, "n_features": numeric_df.shape[1], "n_samples": numeric_df.shape[0]},
    }

    # Optionally save
    if save_path_prefix:
        path = _save_bytes_to_file(png_bytes, save_path_prefix, "png", results_dir=results_dir)
        out["file_path"] = path

    return out


def plot_volcano(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    logfc_col: str = "logFC",
    pval_col: str = "pvalue",
    label_col: Optional[str] = None,
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 1.0,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    """
    Volcano plot. Input df must have logFC and p-value columns.
    """
    return_types = return_types or ["plotly", "png_bytes"]
    df = _read_df(df, csv_path)
    for c in [logfc_col, pval_col]:
        if c not in df.columns:
            raise ValueError(f"Column `{c}` not found in dataframe.")

    # prepare
    eps = np.nextafter(0, 1)
    df = df.copy()
    df["-log10p"] = -np.log10(df[pval_col].replace(0, eps))
    df["significant"] = (df[pval_col] <= pval_cutoff) & (np.abs(df[logfc_col]) >= logfc_cutoff)

    # Plotly
    fig_plotly = None
    if "plotly" in return_types:
        fig_plotly = go.Figure()
        fig_plotly.add_trace(go.Scattergl(
            x=df[logfc_col],
            y=df["-log10p"],
            mode="markers",
            marker=dict(size=6, color="rgba(80,80,80,0.8)"),
            text=df[label_col] if label_col and label_col in df.columns else None,
            name="all"
        ))
        sig = df[df["significant"]]
        if not sig.empty:
            fig_plotly.add_trace(go.Scattergl(
                x=sig[logfc_col],
                y=sig["-log10p"],
                mode="markers",
                marker=dict(size=6, symbol="diamond", color="red"),
                name="significant"
            ))
        fig_plotly.update_layout(title="Volcano plot", xaxis_title="log2(Fold Change)", yaxis_title="-log10(p-value)")

    # Matplotlib static
    fig_mpl = plt.figure(figsize=(6, 4))
    ax = fig_mpl.add_subplot(111)
    ax.scatter(df[logfc_col], df["-log10p"], s=10)
    if not df[df["significant"]].empty:
        ax.scatter(df[df["significant"]][logfc_col], df[df["significant"]]["-log10p"], s=12, marker="D", color="red")
    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Volcano plot")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly) if fig_plotly is not None else None

    sig_labels = []
    if label_col and label_col in df.columns:
        sig_labels = df[df["significant"]][label_col].astype(str).tolist()

    out = {
        "plotly": _fig_to_plotly_json(fig_plotly) if fig_plotly is not None else None,
        "png_bytes": png_bytes,
        "plotly_png_bytes": plotly_png,
        "meta": {"significant_count": int(df["significant"].sum()), "significant_labels": sig_labels},
    }
    if save_path_prefix:
        path = _save_bytes_to_file(png_bytes, save_path_prefix, "png", results_dir=results_dir)
        out["file_path"] = path
    return out


def plot_manhattan(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    pval_col: str = "pvalue",
    snp_col: Optional[str] = "snp",
    significance_line: Optional[float] = None,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    """
    Simple Manhattan plot. Chromosome can be numeric or strings (X, Y).
    """
    return_types = return_types or ["plotly", "png_bytes"]
    df = _read_df(df, csv_path)
    for c in [chrom_col, pos_col, pval_col]:
        if c not in df.columns:
            raise ValueError(f"Column `{c}` not found.")

    dd = df.copy()
    # order chromosomes by natural order (attempt numeric, else string)
    unique_chroms = list(dd[chrom_col].astype(str).unique())
    # map to ordered offsets
    chrom_offsets = {}
    cumulative = 0
    grouped = dd.groupby(chrom_col)
    ordered_chroms = list(grouped.groups.keys())
    for chrom in ordered_chroms:
        grp = grouped.get_group(chrom)
        chrom_offsets[chrom] = cumulative
        cumulative += int(grp[pos_col].max() if grp[pos_col].max() is not None else 0)

    dd["cumulative_pos"] = dd[pos_col] + dd[chrom_col].map(chrom_offsets)
    dd["-log10p"] = -np.log10(dd[pval_col].replace(0, np.nextafter(0, 1)))

    # Plotly
    fig_plotly = None
    if "plotly" in return_types:
        fig_plotly = go.Figure()
        fig_plotly.add_trace(go.Scattergl(
            x=dd["cumulative_pos"],
            y=dd["-log10p"],
            mode="markers",
            marker=dict(size=4),
            text=dd[snp_col] if snp_col and snp_col in dd.columns else None
        ))
        if significance_line:
            fig_plotly.add_hline(y=-np.log10(significance_line), line_dash="dash", line_color="red")
        fig_plotly.update_layout(title="Manhattan plot", xaxis_title="Genomic position", yaxis_title="-log10(p-value)")

    # Matplotlib static
    fig_mpl = plt.figure(figsize=(10, 4))
    ax = fig_mpl.add_subplot(111)
    ax.scatter(dd["cumulative_pos"], dd["-log10p"], s=3)
    if significance_line:
        ax.axhline(-np.log10(significance_line), color="red", linestyle="--")
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Manhattan plot")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly) if fig_plotly is not None else None

    out = {
        "plotly": _fig_to_plotly_json(fig_plotly) if fig_plotly is not None else None,
        "png_bytes": png_bytes,
        "plotly_png_bytes": plotly_png,
        "meta": {"n_snps": int(len(dd))},
    }
    if save_path_prefix:
        path = _save_bytes_to_file(png_bytes, save_path_prefix, "png", results_dir=results_dir)
        out["file_path"] = path
    return out


def plot_heatmap(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    row_names: Optional[str] = None,
    col_names: Optional[str] = None,
    cluster_rows: bool = False,
    cluster_cols: bool = False,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    """
    Simple heatmap. Expects a numeric matrix: rows = features, columns = samples.
    """
    return_types = return_types or ["plotly", "png_bytes"]
    df = _read_df(df, csv_path)

    numeric_df = df.select_dtypes(include=[np.number])
    if numeric_df.shape[0] == 0 or numeric_df.shape[1] == 0:
        raise ValueError("No numeric data available for heatmap.")

    # Plotly heatmap
    fig_plotly = None
    if "plotly" in return_types:
        fig_plotly = go.Figure(data=go.Heatmap(
            z=numeric_df.values,
            x=numeric_df.columns.tolist(),
            y=numeric_df.index.astype(str).tolist(),
            colorbar=dict(title="value")
        ))
        fig_plotly.update_layout(title="Heatmap")

    # Matplotlib static
    fig_mpl = plt.figure(figsize=(8, 6))
    ax = fig_mpl.add_subplot(111)
    cax = ax.imshow(numeric_df.values, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(numeric_df.shape[1]))
    ax.set_xticklabels(numeric_df.columns, rotation=90, fontsize=6)
    ax.set_yticks(np.arange(numeric_df.shape[0]))
    ax.set_yticklabels(numeric_df.index.astype(str), fontsize=6)
    fig_mpl.colorbar(cax)
    ax.set_title("Heatmap")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly) if fig_plotly is not None else None

    out = {
        "plotly": _fig_to_plotly_json(fig_plotly) if fig_plotly is not None else None,
        "png_bytes": png_bytes,
        "plotly_png_bytes": plotly_png,
        "meta": {"shape": numeric_df.shape},
    }
    if save_path_prefix:
        path = _save_bytes_to_file(png_bytes, save_path_prefix, "png", results_dir=results_dir)
        out["file_path"] = path
    return out


def plot_boxplot(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    value_col: Optional[str] = None,
    group_col: Optional[str] = None,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    """
    Boxplot for grouped data or a single distribution.
    - If group_col is provided, draws box per group.
    - If not, draws a single box for value_col.
    """
    return_types = return_types or ["plotly", "png_bytes"]
    df = _read_df(df, csv_path)
    if value_col is None or value_col not in df.columns:
        raise ValueError("`value_col` must be provided and present in the DataFrame.")

    # Plotly
    fig_plotly = None
    if "plotly" in return_types:
        if group_col and group_col in df.columns:
            fig_plotly = px.box(df, x=group_col, y=value_col, title="Boxplot")
        else:
            fig_plotly = px.box(df, y=value_col, title="Boxplot")

    # Matplotlib
    fig_mpl = plt.figure(figsize=(6, 4))
    ax = fig_mpl.add_subplot(111)
    if group_col and group_col in df.columns:
        groups = [g[value_col].dropna().values for _, g in df.groupby(group_col)]
        ax.boxplot(groups)
        ax.set_xticklabels([str(k) for k in df[group_col].unique()], rotation=45)
    else:
        ax.boxplot(df[value_col].dropna().values)
    ax.set_title("Boxplot")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly) if fig_plotly is not None else None

    out = {
        "plotly": _fig_to_plotly_json(fig_plotly) if fig_plotly is not None else None,
        "png_bytes": png_bytes,
        "plotly_png_bytes": plotly_png,
        "meta": {"n_points": int(df[value_col].notna().sum())},
    }
    if save_path_prefix:
        path = _save_bytes_to_file(png_bytes, save_path_prefix, "png", results_dir=results_dir)
        out["file_path"] = path
    return out


def igv_track_meta(name: str, url: str, filetype: Optional[str] = None, index_url: Optional[str] = None) -> Dict[str, Any]:
    """
    Return an IGV.js-compatible track metadata object (this service does not host the file).
    """
    track = {"name": name, "url": url}
    if filetype:
        track["type"] = filetype
    if index_url:
        track["indexURL"] = index_url
    return track


# -------------------------
# Helper: convenience wrapper to return the preferred output
# -------------------------
def choose_best_output(out: Dict[str, Any], prefer: str = "plotly") -> Tuple[str, bytes]:
    """
    Returns (mime_type, bytes) for the best available output.
    prefer: 'plotly' or 'png' or 'plotly_png'
    - 'plotly' returns JSON bytes with mime 'application/json'
    - 'png' returns matplotlib PNG bytes with mime 'image/png'
    - 'plotly_png' attempts to return plotly-exported PNG bytes (kaleido)
    """
    prefer = prefer.lower()
    if prefer == "plotly" and out.get("plotly") is not None:
        return "application/json", json.dumps(out["plotly"]).encode("utf-8")
    if prefer == "plotly_png" and out.get("plotly_png_bytes") is not None:
        return "image/png", out["plotly_png_bytes"]
    # fallback to matplotlib png
    if out.get("png_bytes") is not None:
        return "image/png", out["png_bytes"]
    # last resort: if plotly available, return JSON
    if out.get("plotly") is not None:
        return "application/json", json.dumps(out["plotly"]).encode("utf-8")
    raise RuntimeError("No usable output in provided result.")

# -------------------------
# Additional plotting functions
# -------------------------

def plot_histogram(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    column: str = None,
    bins: int = 30,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    df = _read_df(df, csv_path)
    return_types = return_types or ["plotly", "png_bytes"]
    if column not in df.columns:
        raise ValueError(f"Column `{column}` not found.")

    data = df[column].dropna()
    
    # Plotly
    fig_plotly = px.histogram(data, x=column, nbins=bins, title=f"Histogram of {column}")

    # Matplotlib
    fig_mpl = plt.figure(figsize=(6, 4))
    plt.hist(data, bins=bins, color="skyblue", edgecolor="black")
    plt.xlabel(column)
    plt.ylabel("Count")
    plt.title(f"Histogram of {column}")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly)

    out = {
        "plotly": _fig_to_plotly_json(fig_plotly),
        "png_bytes": png_bytes,
        "plotly_png_bytes": plotly_png,
        "meta": {"n_points": int(len(data))}
    }
    if save_path_prefix:
        out["file_path"] = _save_bytes_to_file(png_bytes, save_path_prefix, results_dir=results_dir)
    return out


def plot_bar(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    category_col: str = None,
    value_col: str = None,
    stacked: bool = False,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    df = _read_df(df, csv_path)
    return_types = return_types or ["plotly", "png_bytes"]
    if category_col not in df.columns or value_col not in df.columns:
        raise ValueError("Columns not found in dataframe.")

    # Prepare grouped data for stacked bar
    if stacked and df[value_col].dtype.kind in 'O':
        df[value_col] = pd.to_numeric(df[value_col], errors='coerce')

    grouped = df.groupby(category_col)[value_col].sum().reset_index()

    # Plotly
    fig_plotly = px.bar(grouped, x=category_col, y=value_col, title=f"Bar plot of {value_col} by {category_col}")

    # Matplotlib
    fig_mpl = plt.figure(figsize=(6, 4))
    plt.bar(grouped[category_col], grouped[value_col], color="skyblue", edgecolor="black")
    plt.xlabel(category_col)
    plt.ylabel(value_col)
    plt.title(f"Bar plot of {value_col} by {category_col}")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly)

    out = {"plotly": _fig_to_plotly_json(fig_plotly),
           "png_bytes": png_bytes,
           "plotly_png_bytes": plotly_png,
           "meta": {"n_categories": len(grouped)}}
    if save_path_prefix:
        out["file_path"] = _save_bytes_to_file(png_bytes, save_path_prefix, results_dir=results_dir)
    return out


def plot_pie_donut(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    column: str = None,
    donut: bool = False,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    df = _read_df(df, csv_path)
    return_types = return_types or ["plotly", "png_bytes"]
    counts = df[column].value_counts()

    # Plotly
    fig_plotly = px.pie(values=counts.values, names=counts.index, title=f"{'Donut' if donut else 'Pie'} chart of {column}")
    if donut:
        fig_plotly.update_traces(hole=0.4)

    # Matplotlib
    fig_mpl = plt.figure(figsize=(6, 6))
    plt.pie(counts.values, labels=counts.index, autopct="%1.1f%%", startangle=90, wedgeprops=dict(width=0.4 if donut else 1))
    plt.title(f"{'Donut' if donut else 'Pie'} chart of {column}")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly)

    out = {"plotly": _fig_to_plotly_json(fig_plotly),
           "png_bytes": png_bytes,
           "plotly_png_bytes": plotly_png,
           "meta": {"n_categories": len(counts)}}
    if save_path_prefix:
        out["file_path"] = _save_bytes_to_file(png_bytes, save_path_prefix, results_dir=results_dir)
    return out


def plot_line_scatter(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    x_col: str = None,
    y_col: str = None,
    kind: str = "scatter",
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    df = _read_df(df, csv_path)
    return_types = return_types or ["plotly", "png_bytes"]
    if x_col not in df.columns or y_col not in df.columns:
        raise ValueError("Columns not found in dataframe.")

    # Plotly
    if kind == "line":
        fig_plotly = px.line(df, x=x_col, y=y_col, title=f"Line plot {y_col} vs {x_col}")
    else:
        fig_plotly = px.scatter(df, x=x_col, y=y_col, title=f"Scatter plot {y_col} vs {x_col}")

    # Matplotlib
    fig_mpl = plt.figure(figsize=(6, 4))
    if kind == "line":
        plt.plot(df[x_col], df[y_col], marker='o', color='skyblue')
    else:
        plt.scatter(df[x_col], df[y_col], color='skyblue')
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f"{kind.capitalize()} plot of {y_col} vs {x_col}")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly)

    out = {"plotly": _fig_to_plotly_json(fig_plotly),
           "png_bytes": png_bytes,
           "plotly_png_bytes": plotly_png,
           "meta": {"n_points": len(df)}}
    if save_path_prefix:
        out["file_path"] = _save_bytes_to_file(png_bytes, save_path_prefix, results_dir=results_dir)
    return out


def plot_lollipop(
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    category_col: str = None,
    value_col: str = None,
    return_types: Optional[List[str]] = None,
    save_path_prefix: Optional[str] = None,
    results_dir: str = DEFAULT_RESULTS_DIR,
) -> Dict[str, Any]:
    """
    Lollipop plot (category vs value)
    """
    df = _read_df(df, csv_path)
    return_types = return_types or ["plotly", "png_bytes"]
    df_sorted = df.sort_values(value_col)

    # Plotly
    fig_plotly = go.Figure()
    fig_plotly.add_trace(go.Scatter(x=df_sorted[category_col], y=df_sorted[value_col], mode='markers+lines', marker=dict(size=10)))
    fig_plotly.update_layout(title=f"Lollipop plot: {value_col} by {category_col}", xaxis_title=category_col, yaxis_title=value_col)

    # Matplotlib
    fig_mpl = plt.figure(figsize=(6, 4))
    plt.stem(df_sorted[category_col], df_sorted[value_col], basefmt=" ", use_line_collection=True)
    plt.xlabel(category_col)
    plt.ylabel(value_col)
    plt.title(f"Lollipop plot: {value_col} by {category_col}")
    png_bytes = _fig_to_png_bytes_matplotlib(fig_mpl)
    plotly_png = _plotly_fig_to_png_bytes(fig_plotly)

    out = {"plotly": _fig_to_plotly_json(fig_plotly),
           "png_bytes": png_bytes,
           "plotly_png_bytes": plotly_png,
           "meta": {"n_categories": len(df_sorted)}}
    if save_path_prefix:
        out["file_path"] = _save_bytes_to_file(png_bytes, save_path_prefix, results_dir=results_dir)
    return out


# -------------------------
# Example usage (for reference)
# -------------------------
if __name__ == "__main__":
    # demo with random data
    import pandas as pd
    rng = np.random.default_rng(1)
    df_demo = pd.DataFrame(rng.normal(size=(50, 10)), columns=[f"g{i}" for i in range(10)])
    df_demo["sample_id"] = [f"s{i}" for i in range(len(df_demo))]

    res = plot_pca(df=df_demo, sample_column="sample_id", save_path_prefix="demo_pca")
    print("Explained variance:", res["explained_variance_ratio"])
    # save output file if present
    if res.get("file_path"):
        print("Saved to:", res["file_path"])
