"""
Demo for BEDParser with CacheManager
"""

import time
from omnibioai.services.annotation_service.cache_utils import CacheManager
from omnibioai.services.annotation_service.parsers.bed_parser import BEDParser

def main():
    bed_file = "data/sample_regions.bed"

    # Initialize BEDParser and CacheManager
    parser = BEDParser(bed_file)
    cache = CacheManager(namespace="bed_demo", use_pickle=True, ttl=10)  # 10-second TTL

    # Single query
    region_key = "chr1:120-180"
    cached_result = cache.get(region_key)
    if cached_result:
        print("Cache hit!")
        print(cached_result)
    else:
        print("Cache miss, queried from BED file:")
        result = parser.query_region("chr1", 120, 180)
        print("Query results for chr1:120-180")
        for r in result:
            print(r)
        cache.set(region_key, result)

    # Batch queries
    batch_regions = [
        ("chr1", 100, 150),
        ("chr1", 140, 210),
        ("chr1", 200, 250)
    ]

    print("\nBatch query results:")
    for i, (chrom, start, end) in enumerate(batch_regions, 1):
        key = f"{chrom}:{start}-{end}"
        cached = cache.get(key)
        if cached:
            print(f"Region {i}: Cache hit!")
        else:
            print(f"Region {i}: Cache miss, queried from BED file")
            cached = parser.query_region(chrom, start, end)
            cache.set(key, cached)
        print(f"Region {i} results: {cached}")

    # TTL expiration demonstration
    print("\nWaiting 12 seconds for TTL expiration...")
    time.sleep(12)

    # Check a previously cached region
    expired_result = cache.get(region_key)
    print(f"After TTL expiration, cached result for {region_key}: {expired_result}")

    # Clear cache
    cache.clear()
    print("Cache cleared. Any query now returns:", cache.get(region_key))


if __name__ == "__main__":
    main()
