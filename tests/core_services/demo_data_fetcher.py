"""
Demo script for DataFetcher

Downloads and caches datasets by specifying name, version, and URL.
"""

from omnibioai.services.dataset_service.data_fetcher import DataFetcher

def main():
    # Initialize fetcher with cache directory
    fetcher = DataFetcher(cache_dir="data/cache")

    # Example: Fetch gnomAD v3.1
    gnomad_v3_1_url = (
        "https://storage.googleapis.com/gnomad-public/release/3.1/vcf/"
        "gnomad.genomes.v3.1.sites.chr1.vcf.bgz"
    )
    dataset_info_v3_1 = fetcher.fetch(name="gnomAD", version="v3.1", url=gnomad_v3_1_url)
    print("gnomAD v3.1:", dataset_info_v3_1)

    # Example: Fetch gnomAD v2.1
    gnomad_v2_1_url = (
        "https://storage.googleapis.com/gnomad-public/release/2.1/vcf/"
        "gnomad.genomes.r2.1.sites.chr1.vcf.gz"
    )
    dataset_info_v2_1 = fetcher.fetch(name="gnomAD", version="v2.1", url=gnomad_v2_1_url)
    print("gnomAD v2.1:", dataset_info_v2_1)

if __name__ == "__main__":
    main()
