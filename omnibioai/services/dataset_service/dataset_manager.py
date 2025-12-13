"""
Module: dataset_manager
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the DatasetManager class for managing datasets in OmniBioAI.
    Supports fetching hosted datasets, local caching, versioning, registration,
    and applying transformations (e.g., normalization, indexing, annotation prep).

Usage:
    from omnibioai.services.dataset_manager import DatasetManager

    manager = DatasetManager()

    # Fetch the latest version of a hosted dataset
    dataset_info = manager.get_dataset("gnomAD")
    print(dataset_info)
    # Output: {'files': ['~/.omnibioai_datasets/gnomAD_v3.1.vcf.gz'], 'version': 'v3.1'}

    # Transform the dataset (e.g., normalization)
    transformed_info = manager.transform_dataset("gnomAD", transformation="normalize")
    print(transformed_info)
    # Output: {'files': ['~/.omnibioai_datasets/gnomAD_v3.1_normalize.vcf.gz'], 'version': 'v3.1'}

Classes:
    - DatasetManager:
        Manages datasets including fetching, caching, versioning, registration,
        and transformations.

        Methods:
            * __init__(cache_dir: Optional[str] = None):
                Initializes the manager with a cache directory, DataFetcher, and DataTransformers.
            * register_dataset(name: str, files: list, version: str):
                Registers a local dataset with name, files, and version.
            * get_dataset(name: str, version: Optional[str] = None) -> dict:
                Retrieves dataset by name and optional version. Downloads if not cached.
            * transform_dataset(name: str, version: Optional[str] = None, transformation: str = "normalize") -> dict:
                Applies a transformation to the dataset and returns transformed file paths.

Constants:
    - DEFAULT_CACHE_DIR: Default cache directory in the user home folder.

Dependencies:
    - os: For path and directory management.
    - typing: For type hints (Optional, Dict, Any).
    - DataFetcher: Downloads and caches hosted datasets.
    - DataTransformers: Applies transformations to datasets.
"""


import os
from typing import Optional, Dict, Any
from .data_fetcher import DataFetcher
from .data_transformers import DataTransformers

DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".omnibioai_datasets")

class DatasetManager:
    """
    Manages datasets: fetching, caching, versioning, and transformations.
    """

    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        os.makedirs(self.cache_dir, exist_ok=True)
        self.fetcher = DataFetcher(self.cache_dir)
        self.transformer = DataTransformers()
        self.datasets: Dict[str, Dict[str, Any]] = {}  # {dataset_name: {version, files}}

    def register_dataset(self, name: str, files: list, version: str):
        """
        Register a dataset locally with a name, files, and version.
        """
        self.datasets[name] = {"files": files, "version": version}

    def get_dataset(self, name: str, version: Optional[str] = None) -> Dict[str, Any]:
        """
        Retrieve dataset by name and optional version.
        Downloads if not cached.
        """
        if name not in self.datasets:
            # Attempt to fetch dataset
            dataset_info = self.fetcher.fetch(name, version)
            self.datasets[name] = dataset_info
        else:
            dataset_info = self.datasets[name]
            if version and dataset_info["version"] != version:
                dataset_info = self.fetcher.fetch(name, version)
                self.datasets[name] = dataset_info

        return dataset_info

    def transform_dataset(self, name: str, version: Optional[str] = None, transformation: str = "normalize"):
        """
        Apply transformation to dataset.
        """
        dataset_info = self.get_dataset(name, version)
        transformed_files = self.transformer.apply(dataset_info["files"], transformation)
        return {"files": transformed_files, "version": dataset_info["version"]}
