# demo_time_series_service.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from omnibioai.services.ml_service.time_series_service import TimeSeriesService
import os
from datetime import datetime

# Ensure reports folder exists
REPORTS_DIR = "data/reports"
os.makedirs(REPORTS_DIR, exist_ok=True)

def save_plot(fig, name: str):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(REPORTS_DIR, f"{name}_{timestamp}.png")
    fig.savefig(path)
    plt.close(fig)
    print(f"[Saved] {path}")

def plot_forecast_with_ci(series, model_result, steps=5, title="Forecast"):
    """
    Plot original series and forecast with 95% confidence interval.
    Works for ARIMAResults and SARIMAXResults.
    """
    fig, ax = plt.subplots(figsize=(10,5))
    series.plot(ax=ax, label="Original")

    # For ARIMA / SARIMA
    if hasattr(model_result, "get_forecast"):
        forecast_obj = model_result.get_forecast(steps=steps)
        forecast = forecast_obj.predicted_mean
        ci = forecast_obj.conf_int()
        ax.plot(forecast.index, forecast, '--', label="Forecast")
        ax.fill_between(forecast.index, ci.iloc[:,0], ci.iloc[:,1], color='orange', alpha=0.2, label="95% CI")
    else:  # Bayesian placeholder
        forecast = model_result.get("forecast", [])
        ax.plot(range(len(series), len(series)+len(forecast)), forecast, '--', label="Forecast")

    ax.set_title(title)
    ax.legend()
    return fig

def main():
    # Sample time series data
    np.random.seed(42)
    data = pd.Series(np.random.randn(24).cumsum() + 100)

    ts_service = TimeSeriesService()

    # ARIMA
    arima_model = ts_service.train_arima(data)
    print("=== ARIMA Forecast ===")
    print(arima_model["forecast"])
    fig = plot_forecast_with_ci(data, arima_model["model"], steps=5, title="ARIMA Forecast")
    save_plot(fig, "arima_forecast")

    # SARIMA
    sarima_model = ts_service.train_sarima(data, seasonal_order=(1,1,1,12))
    print("\n=== SARIMA Forecast ===")
    print(sarima_model["forecast"])
    fig = plot_forecast_with_ci(data, sarima_model["model"], steps=5, title="SARIMA Forecast")
    save_plot(fig, "sarima_forecast")

    # Bayesian placeholder
    bayesian_model = ts_service.train_bayesian(data)
    print("\n=== Bayesian Forecast ===")
    print(bayesian_model["status"])

if __name__ == "__main__":
    main()
