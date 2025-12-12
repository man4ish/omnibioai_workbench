from omnibioai.services.model_zoo_service import ModelZooService
from omnibioai.services.logger_service import logger

# ----------------------------
# 1. Initialize service
# ----------------------------
model_zoo = ModelZooService()

# ----------------------------
# 2. Define test cases
# ----------------------------
test_tasks = {
    "single_cell": ["Scanpy-ML", "scVI", "CellTypist"],
    "variant": ["VEP", "SnpEff", "DeepVariant"],
    "multi_omics": ["MOFA2", "MultiPLIER", "DeepBioML"],
    "unknown_task": ["RandomForest", "XGBoost", "BioBERT"],  # default
}

# ----------------------------
# 3. Run tests
# ----------------------------
for task, expected_models in test_tasks.items():
    models = model_zoo.list_models(task)
    if models == expected_models:
        logger.info(f"Task '{task}' passed: {models}")
    else:
        logger.error(f"Task '{task}' failed: expected {expected_models}, got {models}")

