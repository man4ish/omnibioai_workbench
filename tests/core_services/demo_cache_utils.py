"""
Demo: CacheManager
Author: Manish Kumar
Date: 2025-12-13

Description:
    Demonstrates usage of CacheManager for single entry caching,
    batch caching, retrieval, TTL expiration, and clearing cache.
"""

import time
from omnibioai.services.annotation_service.cache_utils import CacheManager

def main():
    # Initialize cache manager
    cache = CacheManager(namespace="gnomAD_demo", use_pickle=True, ttl=5)  # 5 seconds TTL for demo

    print("=== Single Entry Caching ===")
    cache.set("1:12345:A:T", {"allele_freq": 0.05})
    result = cache.get("1:12345:A:T")
    print("Cached value:", result)

    print("\n=== Batch Caching ===")
    batch_data = {
        "1:12345:A:T": {"allele_freq": 0.05},
        "2:54321:G:C": {"allele_freq": 0.12},
        "3:11111:T:G": {"allele_freq": 0.07}
    }
    cache.set_batch(batch_data)

    keys = ["1:12345:A:T", "2:54321:G:C", "3:11111:T:G", "4:99999:C:A"]
    batch_result = cache.get_batch(keys)
    print("Batch retrieval:", batch_result)

    print("\n=== TTL Expiration Demo ===")
    print("Waiting for 6 seconds to expire cache...")
    time.sleep(6)
    expired_value = cache.get("1:12345:A:T")
    print("After TTL expiration:", expired_value)

    print("\n=== Clear Cache ===")
    cache.clear()
    cleared_value = cache.get("2:54321:G:C")
    print("After clearing cache:", cleared_value)

if __name__ == "__main__":
    main()

