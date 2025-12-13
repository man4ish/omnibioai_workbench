"""
Module: dataset_manager
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the DatasetManager class for managing local reference datasets with versioning.
    Loads a dataset registry from JSON and allows retrieving datasets by name and optional version.
    Defaults to returning the latest version if no version is specified.

Usage:
    from omnibioai.services.dataset_manager import DatasetManager

    manager = DatasetManager(registry_path="data/reference/dataset_registry.json")

    # Get the latest dataset by name
    dataset = manager.get_dataset("gnomAD")

    # Get a specific version
    dataset_v1 = manager.get_dataset("gnomAD", version="v2.1")

Classes:
    - DatasetManager:
        Manages local reference datasets and provides versioned access.

        Methods:
            * __init__(registry_path="data/reference/dataset_registry.json"):
                Initializes the manager by loading the dataset registry JSON.
                Raises FileNotFoundError if the registry file does not exist.
            * get_dataset(name: str, version: str = None):
                Retrieves the dataset by name and optional version.
                Returns the latest dataset if version is None.
                Raises ValueError if dataset or version is not found.

Dependencies:
    - json: For reading the dataset registry file.
    - pathlib.Path: For file path handling.
"""


import json
from pathlib import Path

class DatasetManager:
    """
    Manages local reference datasets with versioning.
    """
    def __init__(self, registry_path="data/reference/dataset_registry.json"):
        self.registry_path = Path(registry_path)
        if not self.registry_path.exists():
            raise FileNotFoundError(f"Dataset registry not found at {self.registry_path}")
        with open(self.registry_path) as f:
            self.registry = json.load(f)

    def get_dataset(self, name: str, version: str = None):
        datasets = self.registry.get(name, [])
        if not datasets:
            raise ValueError(f"Dataset '{name}' not found")
        if version:
            for d in datasets:
                if d["version"] == version:
                    return d
            raise ValueError(f"Version '{version}' of dataset '{name}' not found")
        return datasets[-1]  # return latest by default
