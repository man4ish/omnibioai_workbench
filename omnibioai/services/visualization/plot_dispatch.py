"""
plot_dispatch.py

Centralized plot dispatching service for OmnibioAI.

This module provides a unified interface to generate a wide variety of plots
by dispatching calls to specific plotting functions in `visualization_service`.
It supports both exploratory and publication-ready visualizations.

Supported plot types:
- PCA: Principal Component Analysis plots
- Volcano: Volcano plots for differential expression
- Manhattan: Manhattan plots for GWAS data
- Heatmap: Hierarchical clustering heatmaps
- Boxplot: Box-and-whisker plots
- Histogram: Data distribution histograms
- Scatter: Scatter plots
- Bar: Simple bar plots
- Stacked_bar: Stacked bar plots
- Pie: Pie charts
- Lollipop: Lollipop charts

Functions
---------
plot_dispatch(plot_type, df=None, csv_path=None, save_path_prefix=None,
              return_types=None, results_dir=None, **kwargs)
    Dispatch to the appropriate plotting function based on `plot_type`.
    Returns plotly figures, PNG bytes, and metadata for further use.

Example
-------
>>> from plot_dispatch import plot_dispatch
>>> res = plot_dispatch(plot_type="pca", df=my_dataframe, sample_column="sample_id")
>>> print(res["explained_variance_ratio"])
>>> if res.get("file_path"):
>>>     print("Saved to:", res["file_path"])
"""


from typing import Optional, Dict, Any, Union, List
import pandas as pd

from .visualization_service import (
    plot_pca,
    plot_volcano,
    plot_manhattan,
    plot_heatmap,
    plot_boxplot,
    plot_histogram,
    plot_bar,
    plot_lollipop,
)

# Mapping plot type to function
_PLOT_FN_MAP = {
    "pca": plot_pca,
    "volcano": plot_volcano,
    "manhattan": plot_manhattan,
    "heatmap": plot_heatmap,
    "boxplot": plot_boxplot,
    "histogram": plot_histogram,
    "bar": plot_bar,
    "lollipop": plot_lollipop,
}


def plot_dispatch(
    plot_type: str,
    df: Optional[pd.DataFrame] = None,
    csv_path: Optional[str] = None,
    save_path_prefix: Optional[str] = None,
    return_types: Optional[List[str]] = None,
    results_dir: Optional[str] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Dispatch function to call the correct plot function based on `plot_type`.

    Parameters
    ----------
    plot_type : str
        One of 'pca', 'volcano', 'manhattan', 'heatmap', 'boxplot',
        'histogram', 'scatter', 'bar', 'stacked_bar', 'pie', 'lollipop'.
    df : pd.DataFrame, optional
        Input data for plotting.
    csv_path : str, optional
        Path to CSV file if df is not provided.
    save_path_prefix : str, optional
        Prefix for saving the plot PNG.
    return_types : list of str, optional
        List of output types to return. Defaults to ["plotly", "png_bytes"].
    results_dir : str, optional
        Directory to save plots.
    **kwargs
        Additional parameters for specific plot functions.

    Returns
    -------
    dict
        Plotly figure dict, matplotlib PNG bytes, plotly PNG bytes, and meta info.
    """
    plot_type = plot_type.lower()
    if plot_type not in _PLOT_FN_MAP:
        raise ValueError(f"Unsupported plot_type '{plot_type}'. Supported types: {list(_PLOT_FN_MAP.keys())}")

    plot_fn = _PLOT_FN_MAP[plot_type]

    return plot_fn(
        df=df,
        csv_path=csv_path,
        save_path_prefix=save_path_prefix,
        return_types=return_types,
        results_dir=results_dir,
        **kwargs
    )


# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    import numpy as np

    # Demo dataframe
    rng = np.random.default_rng(42)
    df_demo = pd.DataFrame(rng.normal(size=(50, 10)), columns=[f"g{i}" for i in range(10)])
    df_demo["sample_id"] = [f"s{i}" for i in range(len(df_demo))]

    # PCA example
    res = plot_dispatch(
        plot_type="pca",
        df=df_demo,
        sample_column="sample_id",
        save_path_prefix="demo_pca"
    )
    print("Explained variance:", res.get("explained_variance_ratio"))
    if res.get("file_path"):
        print("Saved to:", res["file_path"])
