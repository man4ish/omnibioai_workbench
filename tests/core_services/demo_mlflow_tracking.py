# tests/core_services/demo_mlflow_tracking.py

from omnibioai.services.ml_service.mlflow_tracking import MLFlowTrackingService

def main():
    mlflow_service = MLFlowTrackingService()

    # Start an experiment with optional tags
    exp_info = mlflow_service.start_experiment(
        name="demo_experiment",
        tags={"project": "omnibioai", "type": "demo"}
    )
    print(f"Experiment started: {exp_info}")

    # Log parameters and metrics
    params = {
        "model_type": "RandomForest",
        "n_estimators": 100,
        "max_depth": 5
    }
    metrics = {
        "accuracy": 0.92,
        "f1_score": 0.90
    }
    mlflow_service.log_params_metrics(params, metrics)
    print("Parameters and metrics logged successfully.")

    # End experiment
    mlflow_service.end_experiment()
    print("Experiment ended.")

if __name__ == "__main__":
    main()
