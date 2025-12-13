"""
Module: data_fetcher
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the DataFetcher class to download and cache hosted datasets with versioning.
    Supports fetching datasets by name and version, automatically downloading them if not
    already cached locally.

Usage:
    from omnibioai.services.data_fetcher import DataFetcher

    fetcher = DataFetcher(cache_dir="data/cache")

    # Fetch the latest gnomAD dataset
    dataset_info = fetcher.fetch("gnomAD")
    print(dataset_info)
    # Output: {"files": ["data/cache/gnomAD_v3.1.vcf.gz"], "version": "v3.1"}

    # Fetch a specific version
    dataset_info_v2 = fetcher.fetch("gnomAD", version="v2.1")
    print(dataset_info_v2)

Classes:
    - DataFetcher:
        Downloads and caches hosted reference datasets.

        Methods:
            * __init__(cache_dir: str):
                Initializes the fetcher and ensures the cache directory exists.
            * fetch(name: str, version: Optional[str] = None) -> dict:
                Returns dataset information including local file paths and version.
                Downloads the dataset if it is not already cached.
            * _download(url: str, dest_path: str):
                Internal method to download a file from a URL to the local cache.

Constants:
    - HOSTED_DATASETS: Placeholder dictionary mapping dataset names and versions to URLs.

Dependencies:
    - os: For path and directory management.
    - typing.Optional, Dict, Any: For type hints.
    - requests: For downloading datasets from URLs.
"""

import os
from typing import Optional, Dict, Any

# Placeholder for hosted dataset URLs
HOSTED_DATASETS = {
    "gnomAD": {
        "v3.1": "https://example.com/gnomad_v3.1.vcf.gz",
        "v2.1": "https://example.com/gnomad_v2.1.vcf.gz"
    },
    "clinvar": {
        "2025-12": "https://example.com/clinvar_2025_12.vcf.gz"
    }
}

class DataFetcher:
    """
    Downloads and caches hosted datasets with versioning.
    """

    def __init__(self, cache_dir: str):
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

    def fetch(self, name: str, version: Optional[str] = None) -> Dict[str, Any]:
        if name not in HOSTED_DATASETS:
            raise ValueError(f"Dataset {name} not found in hosted datasets")

        available_versions = HOSTED_DATASETS[name]
        version = version or max(available_versions.keys())  # use latest
        url = available_versions.get(version)

        # Local file path
        local_file = os.path.join(self.cache_dir, f"{name}_{version}.vcf.gz")

        # Download if file does not exist
        if not os.path.exists(local_file):
            self._download(url, local_file)

        return {"files": [local_file], "version": version}

    def _download(self, url: str, dest_path: str):
        """
        Download file from URL to dest_path
        """
        print(f"Downloading {url} -> {dest_path}")
        # Simple implementation using requests
        import requests
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(dest_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete")
