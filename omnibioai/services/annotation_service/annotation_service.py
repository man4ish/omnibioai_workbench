"""
Module: annotation_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the AnnotationService class for core variant annotation in OmniBioAI.
    Supports querying single variants or batches of variants using persistent caching
    to avoid repeated lookups. Integrates with VCFParser and DatasetManager for
    standardized dataset access and annotation.

Usage:
    from omnibioai.services.annotation_service import AnnotationService
    import pandas as pd

    service = AnnotationService()

    # Annotate a single variant
    annotation = service.annotate_variant("1", 1234567, "A", "T")
    print(annotation)

    # Annotate multiple variants from a DataFrame
    df = pd.DataFrame([
        {"chrom": "1", "pos": 1234567, "ref": "A", "alt": "T"},
        {"chrom": "2", "pos": 7654321, "ref": "G", "alt": "C"}
    ])
    annotated_df = service.annotate_variants_batch(df)
    print(annotated_df)

Classes:
    - AnnotationService:
        Core service for variant annotation with persistent caching.

        Methods:
            * __init__():
                Initializes DatasetManager and persistent CacheManager.
            * annotate_variant(chrom: str, pos: int, ref: str, alt: str, dataset_name="gnomAD", version=None):
                Annotates a single variant using the specified dataset and caches the result.
            * annotate_variants_batch(df: pd.DataFrame, dataset_name="gnomAD", version=None):
                Annotates multiple variants provided in a DataFrame and returns a DataFrame of annotations.

Dependencies:
    - pandas: For batch variant handling.
    - DatasetManager: Provides access to variant annotation datasets.
    - VCFParser: Parses VCF files for variant-level queries.
    - CacheManager: Provides persistent caching for annotation results.
"""

import pandas as pd
from .dataset_manager import DatasetManager
from .parsers.vcf_parser import VCFParser
from .cache_utils import CacheManager

class AnnotationService:
    """
    Core variant annotation service with persistent caching
    """
    def __init__(self):
        self.dataset_manager = DatasetManager()
        self.cache = CacheManager()  # persistent cache

    def annotate_variant(self, chrom: str, pos: int, ref: str, alt: str,
                         dataset_name="gnomAD", version=None):
        dataset = self.dataset_manager.get_dataset(dataset_name, version)
        key = f"{chrom}:{pos}:{ref}:{alt}:{dataset_name}:{dataset['version']}"

        # check persistent cache
        cached = self.cache.get(key, version=dataset['version'])
        if cached:
            return cached

        parser = VCFParser(dataset["files"][0])
        annotation = parser.query_variant(chrom, pos, ref, alt)

        # store in persistent cache
        self.cache.set(key, annotation, version=dataset['version'])
        return annotation

    def annotate_variants_batch(self, df: pd.DataFrame, dataset_name="gnomAD", version=None):
        results = []
        for _, row in df.iterrows():
            ann = self.annotate_variant(
                chrom=row["chrom"],
                pos=row["pos"],
                ref=row["ref"],
                alt=row["alt"],
                dataset_name=dataset_name,
                version=version
            )
            results.append(ann)
        return pd.DataFrame(results)
