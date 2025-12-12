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
