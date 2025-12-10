# omnibioai/services/experiment_tracking_service.py

import os
import json
from .logger_service import logger

class ExperimentTrackingService:
    """
    Minimal Experiment Tracking placeholder.
    Logs suggestions and experiment metadata in JSON files for reproducibility.
    """

    def __init__(self, output_dir="data/experiments"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def log_agentic_ai_suggestions(self, dataset_metadata: dict, suggestion_package: dict):
        """
        Save Agentic AI suggestions for a dataset as JSON.
        """
        dataset_id = dataset_metadata.get("id", "unknown_dataset")
        filename = f"{dataset_id}_agentic_suggestions.json"
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump({
                "dataset_metadata": dataset_metadata,
                "suggestions": suggestion_package
            }, f, indent=2)
        logger.info(f"[ExperimentTracking] Logged Agentic AI suggestions at {path}")
        return path

