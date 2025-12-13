"""
Module: data_fetcher
Author: Manish Kumar
Version: 1.0
Date: 2025-12-13

Description:
    Provides the DataFetcher class to download and cache datasets with versioning.
    All dataset URLs and versions are passed as parameters; nothing is hardcoded.
    The class ensures datasets are downloaded only if not already cached.

Usage:
    from omnibioai.services.dataset_service.data_fetcher import DataFetcher

    fetcher = DataFetcher(cache_dir="data/cache")

    # Fetch a dataset
    dataset_info = fetcher.fetch(
        name="gnomAD",
        version="v3.1",
        url="https://gnomad.broadinstitute.org/downloads"
    )
    print(dataset_info)
    # Output: {"files": ["data/cache/gnomAD_v3.1.vcf.gz"], "version": "v3.1"}

Classes:
    - DataFetcher:
        Downloads and caches datasets with versioning.

        Methods:
            __init__(cache_dir: str):
                Initializes the fetcher and ensures the cache directory exists.

            fetch(name: str, version: str, url: str) -> dict:
                Downloads a dataset from the given URL if not already cached.
                Returns dataset information including local file paths and version.

            _download(url: str, dest_path: str):
                Internal method to download a file from a URL to the local cache.

Dependencies:
    - os: For path and directory management.
    - typing: For type hints (Dict, Any).
    - requests: For downloading datasets from URLs.
"""

import os
from typing import Dict, Any

class DataFetcher:
    """
    Downloads and caches datasets with versioning.
    All dataset URLs and versions are passed as parameters.
    """

    def __init__(self, cache_dir: str):
        """
        Initializes the fetcher and ensures the cache directory exists.

        Args:
            cache_dir (str): Path to local directory for caching datasets.
        """
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

    def fetch(self, name: str, version: str, url: str) -> Dict[str, Any]:
        """
        Fetch dataset by name and version using the provided URL.

        Args:
            name (str): Name of the dataset.
            version (str): Version of the dataset.
            url (str): URL to download the dataset from.

        Returns:
            Dict[str, Any]: Dictionary with local file path(s) and version.
        """
        local_file = os.path.join(self.cache_dir, f"{name}_{version}.vcf.gz")

        if not os.path.exists(local_file):
            self._download(url, local_file)

        return {"files": [local_file], "version": version}

    def _download(self, url: str, dest_path: str):
        """
        Download a file from the URL to the local cache directory.

        Args:
            url (str): URL of the file to download.
            dest_path (str): Local path to save the downloaded file.
        """
        print(f"Downloading {url} -> {dest_path}")
        import requests
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(dest_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete")
