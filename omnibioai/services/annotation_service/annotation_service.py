"""
Module: annotation_service
Author: Manish Kumar
Version: 1.1
Date: 2025-12-13

Description:
    Provides the AnnotationService class for core variant annotation in OmniBioAI.
    Supports querying single variants or batches of variants using persistent caching
    to avoid repeated lookups. Integrates with VCFParser and DatasetManager for
    standardized dataset access and annotation.
"""

import pandas as pd
from .dataset_manager import DatasetManager
from .parsers.vcf_parser import VCFParser
from .cache_utils import CacheManager
from typing import Optional, Dict, Any

class AnnotationService:
    """
    Core variant annotation service with persistent caching
    """

    def __init__(self, cache_ttl: Optional[int] = None, use_pickle: bool = True):
        """
        Initializes DatasetManager and CacheManager.

        Args:
            cache_ttl: Optional TTL (time-to-live in seconds) for cache entries
            use_pickle: Use pickle for storing cached objects
        """
        self.dataset_manager = DatasetManager()
        self.cache = CacheManager(namespace="annotation", use_pickle=use_pickle, ttl=cache_ttl)

    def annotate_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        dataset_name: str = "gnomAD",
        version: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Annotates a single variant and caches the result.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            dataset_name: Dataset to use for annotation
            version: Dataset version

        Returns:
            Annotation dictionary
        """
        dataset = self.dataset_manager.get_dataset(dataset_name, version)
        key = f"{chrom}:{pos}:{ref}:{alt}:{dataset_name}:{dataset['version']}"

        cached = self.cache.get(key)
        if cached:
            return cached

        parser = VCFParser(dataset["files"][0])
        annotation = parser.query_variant(chrom, pos, ref, alt)

        self.cache.set(key, annotation)
        return annotation

    def annotate_variants_batch(
        self,
        df: pd.DataFrame,
        dataset_name: str = "gnomAD",
        version: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Annotates multiple variants in a DataFrame efficiently.

        Args:
            df: DataFrame with columns ["chrom", "pos", "ref", "alt"]
            dataset_name: Dataset to use
            version: Dataset version

        Returns:
            DataFrame with annotations
        """
        dataset = self.dataset_manager.get_dataset(dataset_name, version)
        parser = VCFParser(dataset["files"][0])

        results = []
        for _, row in df.iterrows():
            key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}:{dataset_name}:{dataset['version']}"
            cached = self.cache.get(key)
            if cached:
                results.append(cached)
            else:
                annotation = parser.query_variant(row["chrom"], row["pos"], row["ref"], row["alt"])
                self.cache.set(key, annotation)
                results.append(annotation)

        return pd.DataFrame(results)


    """
Module: annotation_service
Author: Manish Kumar
Version: 1.1
Date: 2025-12-13

Description:
    Provides the AnnotationService class for core variant annotation in OmniBioAI.
    Supports querying single variants or batches of variants using persistent caching
    to avoid repeated lookups. Integrates with VCFParser and DatasetManager for
    standardized dataset access and annotation.
"""

import pandas as pd
from .dataset_manager import DatasetManager
from .parsers.vcf_parser import VCFParser
from .cache_utils import CacheManager
from typing import Optional, Dict, Any

class AnnotationService:
    """
    Core variant annotation service with persistent caching
    """

    def __init__(self, cache_ttl: Optional[int] = None, use_pickle: bool = True):
        """
        Initializes DatasetManager and CacheManager.

        Args:
            cache_ttl: Optional TTL (time-to-live in seconds) for cache entries
            use_pickle: Use pickle for storing cached objects
        """
        self.dataset_manager = DatasetManager()
        self.cache = CacheManager(namespace="annotation", use_pickle=use_pickle, ttl=cache_ttl)

    def annotate_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        dataset_name: str = "gnomAD",
        version: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Annotates a single variant and caches the result.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            dataset_name: Dataset to use for annotation
            version: Dataset version

        Returns:
            Annotation dictionary
        """
        dataset = self.dataset_manager.get_dataset(dataset_name, version)
        key = f"{chrom}:{pos}:{ref}:{alt}:{dataset_name}:{dataset['version']}"

        cached = self.cache.get(key)
        if cached:
            return cached

        parser = VCFParser(dataset["files"][0])
        annotation = parser.query_variant(chrom, pos, ref, alt)

        self.cache.set(key, annotation)
        return annotation

    def annotate_variants_batch(
        self,
        df: pd.DataFrame,
        dataset_name: str = "gnomAD",
        version: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Annotates multiple variants in a DataFrame efficiently.

        Args:
            df: DataFrame with columns ["chrom", "pos", "ref", "alt"]
            dataset_name: Dataset to use
            version: Dataset version

        Returns:
            DataFrame with annotations
        """
        dataset = self.dataset_manager.get_dataset(dataset_name, version)
        parser = VCFParser(dataset["files"][0])

        results = []
        for _, row in df.iterrows():
            key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}:{dataset_name}:{dataset['version']}"
            cached = self.cache.get(key)
            if cached:
                results.append(cached)
            else:
                annotation = parser.query_variant(row["chrom"], row["pos"], row["ref"], row["alt"])
                self.cache.set(key, annotation)
                results.append(annotation)

        return pd.DataFrame(results)
    
