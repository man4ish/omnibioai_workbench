"""
Module: cache_utils
Author: Manish Kumar
Version: 2.0
Date: 2025-12-13

Description:
    Provides CacheManager class for advanced persistent caching of annotation results.
    Features:
        - Batch caching
        - Pickle or JSON storage
        - Optional TTL (expiration)
        - Namespace support for separating caches

Usage:
    from omnibioai.services.annotation_service.cache_utils import CacheManager

    cache = CacheManager(namespace="gnomAD", use_pickle=True, ttl=3600)

    # Set single cache entry
    cache.set("1:12345:A:T", {"allele_freq": 0.05})

    # Get cached value
    result = cache.get("1:12345:A:T")
    print(result)

    # Batch cache multiple entries
    batch_data = {
        "1:12345:A:T": {"allele_freq": 0.05},
        "2:54321:G:C": {"allele_freq": 0.12}
    }
    cache.set_batch(batch_data)

    # Clear all cache in namespace
    cache.clear()
"""

import os
import json
import pickle
import time
from typing import Any, Optional, Dict

DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".omnibioai_cache")

class CacheManager:
    """
    Advanced persistent cache manager with versioning, batch caching,
    TTL, and pickle support.
    """

    def __init__(self, namespace: Optional[str] = None, cache_dir: Optional[str] = None,
                 use_pickle: bool = False, ttl: Optional[int] = None):
        """
        Initialize CacheManager.

        Args:
            namespace: Optional namespace to separate cache files (e.g., dataset name)
            cache_dir: Directory for cache storage
            use_pickle: If True, store cache as pickle instead of JSON
            ttl: Time-to-live in seconds; cached entries older than TTL are invalid
        """
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        self.namespace = namespace or "default"
        self.use_pickle = use_pickle
        self.ttl = ttl  # seconds
        self.ns_dir = os.path.join(self.cache_dir, self.namespace)
        os.makedirs(self.ns_dir, exist_ok=True)

    def _get_cache_file(self, key: str) -> str:
        safe_key = key.replace(":", "_").replace("/", "_")
        ext = ".pkl" if self.use_pickle else ".json"
        return os.path.join(self.ns_dir, f"{safe_key}{ext}")

    def _is_expired(self, file_path: str) -> bool:
        if self.ttl is None:
            return False
        age = time.time() - os.path.getmtime(file_path)
        return age > self.ttl

    def set(self, key: str, value: Any):
        """
        Store a single cache entry.
        """
        file_path = self._get_cache_file(key)
        if self.use_pickle:
            with open(file_path, "wb") as f:
                pickle.dump(value, f)
        else:
            with open(file_path, "w") as f:
                json.dump(value, f)

    def get(self, key: str) -> Optional[Any]:
        """
        Retrieve a cache entry. Returns None if not found or expired.
        """
        file_path = self._get_cache_file(key)
        if not os.path.exists(file_path) or self._is_expired(file_path):
            return None
        if self.use_pickle:
            with open(file_path, "rb") as f:
                return pickle.load(f)
        else:
            with open(file_path, "r") as f:
                return json.load(f)

    def set_batch(self, data: Dict[str, Any]):
        """
        Store multiple cache entries at once.
        """
        for key, value in data.items():
            self.set(key, value)

    def get_batch(self, keys: list) -> Dict[str, Any]:
        """
        Retrieve multiple cache entries. Returns dict with missing keys as None.
        """
        return {key: self.get(key) for key in keys}

    def clear(self):
        """
        Clear all cache entries in this namespace.
        """
        for f in os.listdir(self.ns_dir):
            file_path = os.path.join(self.ns_dir, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
