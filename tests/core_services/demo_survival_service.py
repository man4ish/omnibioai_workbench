# tests/core_services/demo_survival_service.py

import os
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from omnibioai.services.ml_service.survival_service import SurvivalService

# Create reports folder if not exists
REPORTS_DIR = "data/reports"
os.makedirs(REPORTS_DIR, exist_ok=True)

def save_plot(fig, name: str):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(REPORTS_DIR, f"{name}_{timestamp}.png")
    fig.savefig(filename)
    plt.close(fig)
    print(f"[Saved] {filename}")

def main():
    service = SurvivalService()

    # -----------------------------
    # Kaplan-Meier example
    # -----------------------------
    durations = pd.Series([5, 6, 6, 2, 4, 4, 8, 9, 3, 10])
    events = pd.Series([1, 0, 0, 1, 1, 1, 0, 0, 1, 0])

    km_result = service.fit_kaplan_meier(durations, events)
    kmf = km_result["model"]

    print("=== Kaplan-Meier Example ===")
    print("Status:", km_result["status"])

    # Plot Kaplan-Meier curve
    fig, ax = plt.subplots()
    kmf.plot_survival_function(ax=ax)
    ax.set_title("Kaplan-Meier Survival Curve")
    save_plot(fig, "kaplan_meier_curve")

    # -----------------------------
    # Cox Regression example
    # -----------------------------
    df = pd.DataFrame({
        "age": [50, 60, 65, 70, 80, 55, 45, 50, 60, 75],
        "treatment": [0, 1, 0, 1, 1, 0, 0, 1, 0, 1],
        "duration": durations,
        "event": events
    })

    cox_result = service.fit_cox_regression(df, duration_col="duration", event_col="event")
    cph = cox_result["model"]

    print("\n=== Cox Regression Example ===")
    print("Status:", cox_result["status"])
    print(cph.summary)

    # Plot partial effects for all numeric covariates
    numeric_covariates = df.select_dtypes(include="number").columns.drop(["duration", "event"])
    fig, ax = plt.subplots(figsize=(8,6))
    for covariate in numeric_covariates:
        values = df[covariate].unique()
        cph.plot_partial_effects_on_outcome(covariates=covariate, values=values, ax=ax, cmap="tab10")
    ax.set_title("Cox Regression: Partial Effects of Numeric Covariates")
    save_plot(fig, "cox_regression_numeric_effects")


if __name__ == "__main__":
    main()
