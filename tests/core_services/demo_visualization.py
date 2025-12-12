# tests/core_services/demo_visualization.py

import os
import numpy as np
import pandas as pd

from omnibioai.services.visualization_service import (
    plot_pca, plot_volcano, plot_manhattan, plot_heatmap, plot_boxplot
)
from omnibioai.services.visualization_ai_client import VisualizationAIClient

# -------------------------
# Setup
# -------------------------
RESULTS_DIR = "tests/core_services/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

rng = np.random.default_rng(42)

# -------------------------
# 1. PCA demo
# -------------------------
df_pca = pd.DataFrame(rng.normal(size=(50, 10)), columns=[f"gene{i}" for i in range(10)])
df_pca["sample_id"] = [f"s{i}" for i in range(len(df_pca))]

print("\n--- PCA Demo ---")
pca_out = plot_pca(df=df_pca, sample_column="sample_id", save_path_prefix="demo_pca", results_dir=RESULTS_DIR)
print("Explained variance:", pca_out.get("explained_variance_ratio"))
if "file_path" in pca_out:
    print("Saved PCA PNG to:", pca_out["file_path"])

# -------------------------
# 2. Boxplot demo
# -------------------------
df_box = pd.DataFrame({
    "value": rng.normal(size=20),
    "group": ["A", "B"] * 10
})

print("\n--- Boxplot Demo ---")
box_out = plot_boxplot(df=df_box, value_col="value", group_col="group",
                       save_path_prefix="demo_boxplot", results_dir=RESULTS_DIR)
if "file_path" in box_out:
    print("Saved Boxplot PNG to:", box_out["file_path"])

# -------------------------
# 3. Heatmap demo
# -------------------------
df_heat = pd.DataFrame(rng.normal(size=(10, 5)), index=[f"gene{i}" for i in range(10)])
print("\n--- Heatmap Demo ---")
heat_out = plot_heatmap(df=df_heat, save_path_prefix="demo_heatmap", results_dir=RESULTS_DIR)
if "file_path" in heat_out:
    print("Saved Heatmap PNG to:", heat_out["file_path"])

# -------------------------
# 4. Volcano plot demo
# -------------------------
df_volcano = pd.DataFrame({
    "logFC": rng.normal(size=50),
    "pvalue": rng.uniform(0, 1, size=50),
    "gene": [f"gene{i}" for i in range(50)]
})

print("\n--- Volcano Plot Demo ---")
volcano_out = plot_volcano(df=df_volcano, logfc_col="logFC", pval_col="pvalue",
                           label_col="gene", save_path_prefix="demo_volcano", results_dir=RESULTS_DIR)
if "file_path" in volcano_out:
    print("Saved Volcano PNG to:", volcano_out["file_path"])

# -------------------------
# 5. AI-assisted visualization demo
# -------------------------
print("\n--- AI-assisted Visualization Demo ---")
viz_ai = VisualizationAIClient(results_dir=RESULTS_DIR)

# AI suggests plot automatically (goal="variance" -> PCA)
ai_out = viz_ai.create_plot(df_pca, goal="variance", save_path_prefix="ai_pca_demo")
print("AI suggested PCA explained variance:", ai_out.get("explained_variance_ratio"))
if "file_path" in ai_out:
    print("Saved AI PCA PNG to:", ai_out["file_path"])

# AI explicitly creates heatmap
ai_heat = viz_ai.create_plot(df_heat, plot_type="heatmap", save_path_prefix="ai_heatmap_demo")
if "file_path" in ai_heat:
    print("Saved AI Heatmap PNG to:", ai_heat["file_path"])

# -------------------------
# Optional: Manhattan demo
# -------------------------
df_manhattan = pd.DataFrame({
    "chrom": ["1"]*10 + ["2"]*10,
    "pos": list(range(1, 21)),
    "pvalue": rng.uniform(0, 1, size=20),
    "snp": [f"rs{i}" for i in range(20)]
})
print("\n--- Manhattan Plot Demo ---")
manhattan_out = plot_manhattan(df=df_manhattan, save_path_prefix="demo_manhattan", results_dir=RESULTS_DIR)
if "file_path" in manhattan_out:
    print("Saved Manhattan PNG to:", manhattan_out["file_path"])
