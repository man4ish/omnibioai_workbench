"""
Demo for OmnibioAI plot_dispatch.py

This script demonstrates the usage of the centralized plotting dispatch service.
Outputs include Plotly figures, PNGs, and metadata saved to the default results directory.
"""

import os
import pandas as pd
import numpy as np
from omnibioai.services.visualization.plot_dispatch import plot_dispatch

# -------------------------
# Ensure results directory
# -------------------------
RESULTS_DIR = "tests/coreservices/demo_outputs"
os.makedirs(RESULTS_DIR, exist_ok=True)

# -------------------------
# Generate sample data
# -------------------------
rng = np.random.default_rng(42)

# PCA / general numeric
df_numeric = pd.DataFrame(rng.normal(size=(50, 5)), columns=[f"gene_{i}" for i in range(5)])
df_numeric["sample_id"] = [f"s{i}" for i in range(len(df_numeric))]

# Volcano plot
df_volcano = pd.DataFrame({
    "gene": [f"g{i}" for i in range(50)],
    "logFC": rng.normal(size=50),
    "pvalue": rng.uniform(0, 1, size=50)
})

# Manhattan plot
df_manhattan = pd.DataFrame({
    "chrom": [f"chr{i%5+1}" for i in range(50)],
    "pos": np.arange(1, 51),
    "pvalue": rng.uniform(0, 1, size=50),
    "snp": [f"rs{i}" for i in range(50)]
})

# Heatmap
df_heatmap = pd.DataFrame(rng.normal(size=(10, 10)), index=[f"gene{i}" for i in range(10)], columns=[f"s{i}" for i in range(10)])

# Boxplot
df_box = pd.DataFrame({
    "group": np.random.choice(["A", "B", "C"], size=50),
    "value": rng.normal(size=50)
})

# Histogram
df_hist = pd.DataFrame({"value": rng.normal(size=100)})

# Bar plot
df_bar = pd.DataFrame({
    "category": np.random.choice(["X", "Y", "Z"], size=50),
    "value": rng.integers(1, 10, size=50)
})

# Pie / Donut
df_pie = pd.DataFrame({"category": np.random.choice(["Cat1", "Cat2", "Cat3"], size=50)})

# Line / Scatter
df_line = pd.DataFrame({
    "x": np.arange(20),
    "y": rng.normal(size=20)
})

# Lollipop
df_lolli = pd.DataFrame({
    "category": [f"C{i}" for i in range(10)],
    "value": rng.normal(size=10)
})

# -------------------------
# Run demo plots
# -------------------------
def run_demo():
    print("=== Running OmnibioAI Plot Dispatch Demo ===")

    # PCA
    res_pca = plot_dispatch(
        plot_type="pca",
        df=df_numeric,
        sample_column="sample_id",
        save_path_prefix="demo_pca/demo_pca",
        results_dir=RESULTS_DIR
    )
    print("PCA explained variance:", res_pca.get("explained_variance_ratio"))

    # Volcano
    res_volcano = plot_dispatch(
        plot_type="volcano",
        df=df_volcano,
        logfc_col="logFC",
        pval_col="pvalue",
        label_col="gene",
        save_path_prefix="demo_volcano/demo_volcano",
        results_dir=RESULTS_DIR
    )
    print("Volcano significant count:", res_volcano["meta"]["significant_count"])

    # Manhattan
    res_manhattan = plot_dispatch(
        plot_type="manhattan",
        df=df_manhattan,
        chrom_col="chrom",
        pos_col="pos",
        pval_col="pvalue",
        snp_col="snp",
        save_path_prefix="demo_manhattan/demo_manhattan",
        results_dir=RESULTS_DIR
    )
    print("Manhattan SNPs:", res_manhattan["meta"]["n_snps"])

    # Heatmap
    res_heatmap = plot_dispatch(
        plot_type="heatmap",
        df=df_heatmap,
        save_path_prefix="demo_heatmap/demo_heatmap",
        results_dir=RESULTS_DIR
    )
    print("Heatmap shape:", res_heatmap["meta"]["shape"])

    # Boxplot
    res_box = plot_dispatch(
        plot_type="boxplot",
        df=df_box,
        value_col="value",
        group_col="group",
        save_path_prefix="demo_boxplot/demo_boxplot",
        results_dir=RESULTS_DIR
    )
    print("Boxplot n_points:", res_box["meta"]["n_points"])

    # Histogram
    res_hist = plot_dispatch(
        plot_type="histogram",
        df=df_hist,
        column="value",
        save_path_prefix="demo_histogram/demo_histogram",
        results_dir=RESULTS_DIR
    )
    print("Histogram n_points:", res_hist["meta"]["n_points"])

    # Bar plot
    res_bar = plot_dispatch(
        plot_type="bar",
        df=df_bar,
        category_col="category",
        value_col="value",
        save_path_prefix="demo_bar/demo_bar",
        results_dir=RESULTS_DIR
    )
    print("Bar n_categories:", res_bar["meta"]["n_categories"])

    # Lollipop
    res_lolli = plot_dispatch(
        plot_type="lollipop",
        df=df_lolli,
        category_col="category",
        value_col="value",
        save_path_prefix="demo_lollipop/demo_lollipop",
        results_dir=RESULTS_DIR
    )
    print("Lollipop n_categories:", res_lolli["meta"]["n_categories"])

    print("=== Demo Finished ===")
    print(f"All outputs saved in: {RESULTS_DIR}")


if __name__ == "__main__":
    run_demo()
