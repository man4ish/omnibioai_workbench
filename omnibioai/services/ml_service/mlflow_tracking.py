# ml_service/mlflow_tracking.py
"""
mlflow_tracking.py

Utility service for integrating MLflow experiment tracking into machine learning workflows.

This module provides a simple wrapper around MLflow to:

- Start and manage experiments
- Log parameters and metrics
- End experiments

Classes
-------
MLFlowTrackingService
    Wrapper for MLflow functionality to track experiments and log results.

Usage Example
-------------
from ml_service import mlflow_tracking

mlflow_service = mlflow_tracking.MLFlowTrackingService()

# Start an experiment
exp_info = mlflow_service.start_experiment("my_experiment")

# Log parameters and metrics
params = {"model_type": "RandomForest", "n_estimators": 100}
metrics = {"accuracy": 0.92, "f1_score": 0.90}
mlflow_service.log_params_metrics(params, metrics)

# End the experiment
mlflow_service.end_experiment()
"""

import mlflow
from typing import Dict, Any

class MLFlowTrackingService:
    """
    MLflow integration for tracking experiments.
    """

    def start_experiment(self, name: str) -> Dict[str, Any]:
        mlflow.set_experiment(name)
        run = mlflow.start_run()
        return {"run_id": run.info.run_id, "status": "started"}

    def log_params_metrics(self, params: Dict[str, Any], metrics: Dict[str, float]) -> None:
        for k, v in params.items():
            mlflow.log_param(k, v)
        for k, v in metrics.items():
            mlflow.log_metric(k, v)

    def end_experiment(self) -> None:
        mlflow.end_run()
