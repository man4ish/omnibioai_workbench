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

# Start an experiment (optional tags)
exp_info = mlflow_service.start_experiment("my_experiment", tags={"project": "omnibioai"})

# Log parameters and metrics
params = {"model_type": "RandomForest", "n_estimators": 100}
metrics = {"accuracy": 0.92, "f1_score": 0.90}
mlflow_service.log_params_metrics(params, metrics)

# End the experiment
mlflow_service.end_experiment()
"""

import mlflow
from typing import Dict, Any, Optional
from datetime import datetime

class MLFlowTrackingService:
    """
    MLflow integration for tracking experiments with optional timestamped run names.
    """

    def start_experiment(self, name: str, tags: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
        """
        Start an MLflow experiment. Adds a timestamp to run name to ensure uniqueness.
        Args:
            name (str): Experiment name.
            tags (dict, optional): Optional tags for the run.
        Returns:
            dict: Run information (run_id and status)
        """
        mlflow.set_experiment(name)
        run_name = f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        run = mlflow.start_run(run_name=run_name, tags=tags)
        return {"run_id": run.info.run_id, "status": "started", "run_name": run_name}

    def log_params_metrics(self, params: Dict[str, Any], metrics: Dict[str, float]) -> None:
        """
        Log parameters and metrics to the current MLflow run.
        """
        for k, v in params.items():
            mlflow.log_param(k, v)
        for k, v in metrics.items():
            mlflow.log_metric(k, v)

    def end_experiment(self) -> None:
        """
        End the current MLflow run.
        """
        mlflow.end_run()
