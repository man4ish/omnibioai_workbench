"""
Module: dataset_manager
Author: Manish Kumar
Version: 1.1
Date: 2025-12-13

Description:
    Provides the DatasetManager class for managing local reference datasets with versioning.
    Loads a dataset registry from JSON and allows retrieving datasets by name and optional version.
    Defaults to returning the latest version if no version is specified. Validates registry format
    and caches lookups for efficiency.

Usage:
    from omnibioai.services.dataset_manager import DatasetManager

    manager = DatasetManager(registry_path="data/reference/dataset_registry.json")

    # Get the latest dataset by name
    dataset = manager.get_dataset("gnomAD")

    # Get a specific version
    dataset_v1 = manager.get_dataset("gnomAD", version="v2.1")
"""

import json
from pathlib import Path
from typing import Optional, Dict, Any
import copy

class DatasetManager:
    """
    Manages local reference datasets with versioning, validation, and caching.
    """
    def __init__(self, registry_path: str = "data/reference/dataset_registry.json"):
        self.registry_path = Path(registry_path)
        if not self.registry_path.exists():
            raise FileNotFoundError(f"Dataset registry not found at {self.registry_path}")

        with open(self.registry_path) as f:
            self.registry: Dict[str, list[Dict[str, Any]]] = json.load(f)

        # Validate registry format
        for name, datasets in self.registry.items():
            if not isinstance(datasets, list):
                raise ValueError(f"Registry entry for '{name}' should be a list")
            for d in datasets:
                if "version" not in d or "files" not in d:
                    raise ValueError(f"Dataset '{name}' entry missing 'version' or 'files': {d}")

        # Simple in-memory cache for faster repeated lookups
        self._cache: Dict[tuple[str, Optional[str]], dict] = {}

    def get_dataset(self, name: str, version: Optional[str] = None) -> dict:
        """
        Retrieve dataset by name and optional version.
        Returns latest version if `version` is None.
        Raises ValueError if dataset or version is not found.
        """
        cache_key = (name, version)
        if cache_key in self._cache:
            return copy.deepcopy(self._cache[cache_key])

        datasets = self.registry.get(name, [])
        if not datasets:
            raise ValueError(f"Dataset '{name}' not found")

        # Retrieve specific version if provided
        if version:
            for d in datasets:
                if d["version"] == version:
                    self._cache[cache_key] = d
                    return copy.deepcopy(d)
            raise ValueError(f"Version '{version}' of dataset '{name}' not found")

        # Return the latest dataset by default (assumes last entry is latest)
        latest = datasets[-1]
        self._cache[cache_key] = latest
        return copy.deepcopy(latest)
