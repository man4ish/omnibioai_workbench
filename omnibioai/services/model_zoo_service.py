"""
Module: model_zoo_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides a minimal Model Zoo service for OmniBioAI. 
    Returns pre-defined recommended models based on task type.
    Useful for suggesting appropriate analysis or ML models for 
    different bioinformatics workflows.

Usage:
    from omnibioai.services.model_zoo_service import ModelZooService

    zoo = ModelZooService()
    models = zoo.list_models("single_cell")
    print(models)  # Output: ["Scanpy-ML", "scVI", "CellTypist"]

Classes:
    - ModelZooService:
        Minimal service to retrieve recommended models for different bioinformatics tasks.
        
        Methods:
            * __init__():
                Initializes the service with a pre-defined mapping of task types to models.
            * list_models(task_type: str) -> List[str]:
                Returns recommended models for the given task type. Defaults to 'general' if unknown.

Dependencies:
    - omnibioai.services.logger_service.logger: For logging recommended models.
"""


from .logger_service import logger

class ModelZooService:
    """
    Minimal Model Zoo placeholder.
    Returns pre-defined models based on task type.
    """

    def __init__(self):
        # Example: pre-defined mapping of tasks to models
        self.models = {
            "single_cell": ["Scanpy-ML", "scVI", "CellTypist"],
            "variant": ["VEP", "SnpEff", "DeepVariant"],
            "multi_omics": ["MOFA2", "MultiPLIER", "DeepBioML"],
            "general": ["RandomForest", "XGBoost", "BioBERT"]
        }

    def list_models(self, task_type: str):
        task_type = task_type.lower()
        recommended = self.models.get(task_type, self.models["general"])
        logger.info(f"[ModelZoo] Recommended models for task '{task_type}': {recommended}")
        return recommended

