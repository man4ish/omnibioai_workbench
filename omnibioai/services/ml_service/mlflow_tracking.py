# ml_service/mlflow_tracking.py
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
