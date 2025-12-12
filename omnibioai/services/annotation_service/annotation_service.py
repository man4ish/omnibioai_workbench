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
