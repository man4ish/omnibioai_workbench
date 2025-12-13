"""
Demo for OmnibioAI visualization_ai_client.py

This script demonstrates the usage of VisualizationAIClient.
It shows how to get AI-assisted plot suggestions and generate plots.
"""

import os
import pandas as pd
import numpy as np

from omnibioai.services.visualization.visualization_ai_client import VisualizationAIClient

# -------------------------
# Setup results directory
# -------------------------
RESULTS_DIR = "tests/coreservices/demo_outputs_ai"
os.makedirs(RESULTS_DIR, exist_ok=True)

# -------------------------
# Generate sample data
# -------------------------
rng = np.random.default_rng(42)

# Numeric + sample ID for PCA/variance
df_numeric = pd.DataFrame(rng.normal(size=(50, 5)), columns=[f"gene_{i}" for i in range(5)])
df_numeric["sample_id"] = [f"s{i}" for i in range(len(df_numeric))]

# Numeric + categorical for boxplot
df_box = pd.DataFrame({
    "group": np.random.choice(["A", "B", "C"], size=50),
    "value": rng.normal(size=50)
})

# Single numeric column for histogram
df_hist = pd.DataFrame({"value": rng.normal(size=100)})

# -------------------------
# Run AI-assisted visualization demo
# -------------------------
def run_demo():
    print("=== Running OmnibioAI AI Visualization Demo ===")
    client = VisualizationAIClient(results_dir=RESULTS_DIR)

    # 1️⃣ PCA suggestion (variance/dimensionality goal)
    plot_type, plot_params = client.suggest_plot(df_numeric, goal="variance")
    print("Suggested plot type for variance goal:", plot_type)
    res_pca = client.create_plot(df_numeric, plot_type=plot_type, save_path_prefix="demo_ai_pca", **plot_params)
    print("PCA explained variance:", res_pca.get("explained_variance_ratio"))
    print("Saved to:", res_pca.get("file_path"))

    # 2️⃣ Boxplot suggestion
    plot_type, plot_params = client.suggest_plot(df_box)
    print("Suggested plot type for boxplot dataset:", plot_type)
    res_box = client.create_plot(df_box, plot_type=plot_type, save_path_prefix="demo_ai_boxplot", **plot_params)
    print("Boxplot n_points:", res_box["meta"]["n_points"])
    print("Saved to:", res_box.get("file_path"))

    # 3️⃣ Histogram suggestion
    plot_type, plot_params = client.suggest_plot(df_hist)
    print("Suggested plot type for histogram dataset:", plot_type)
    res_hist = client.create_plot(df_hist, plot_type=plot_type, save_path_prefix="demo_ai_hist", **plot_params)
    print("Histogram n_points:", res_hist["meta"]["n_points"])
    print("Saved to:", res_hist.get("file_path"))

    print("=== Demo Finished ===")
    print(f"All outputs saved in: {RESULTS_DIR}")


if __name__ == "__main__":
    run_demo()

