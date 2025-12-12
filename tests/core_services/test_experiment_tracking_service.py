import os
import json
from omnibioai.services.experiment_tracking_service import ExperimentTrackingService
from omnibioai.services.logger_service import logger

# ----------------------------
# 1. Initialize service
# ----------------------------
test_output_dir = "data/test_experiments"
exp_service = ExperimentTrackingService(output_dir=test_output_dir)

# ----------------------------
# 2. Prepare test data
# ----------------------------
dataset_metadata = {"id": "dataset_001", "name": "Test Dataset"}
suggestion_package = {"steps": ["normalize", "filter_outliers", "impute_missing"]}

# ----------------------------
# 3. Log suggestions
# ----------------------------
logged_file_path = exp_service.log_agentic_ai_suggestions(dataset_metadata, suggestion_package)
logger.info(f"Logged file path: {logged_file_path}")

# ----------------------------
# 4. Verify file exists and content
# ----------------------------
if os.path.exists(logged_file_path):
    logger.info(f"File successfully created: {logged_file_path}")
    with open(logged_file_path, "r") as f:
        data = json.load(f)
        assert data["dataset_metadata"] == dataset_metadata, "Dataset metadata mismatch"
        assert data["suggestions"] == suggestion_package, "Suggestions mismatch"
        logger.info("JSON content verified successfully")
else:
    logger.error("Failed to create log file")

# ----------------------------
# 5. Cleanup (optional)
# ----------------------------
import shutil
shutil.rmtree(test_output_dir)
logger.info(f"Cleaned up test directory: {test_output_dir}")

