"""
Demo: DatasetService usage
- Fetch hosted datasets (e.g., gnomAD)
- Apply transformations (e.g., normalization)
- Print dataset paths and versions
"""

from omnibioai.core.dataset_service import DatasetManager

def main():
    # Initialize the dataset manager
    dm = DatasetManager()

    # 1. Fetch the latest gnomAD dataset
    print("Fetching gnomAD dataset...")
    gnomad_data = dm.get_dataset("gnomAD")
    print(f"gnomAD version: {gnomad_data['version']}")
    print("Files:", gnomad_data["files"])

    # 2. Fetch a specific version of ClinVar
    print("\nFetching ClinVar dataset version 2025-12...")
    clinvar_data = dm.get_dataset("clinvar", version="2025-12")
    print(f"ClinVar version: {clinvar_data['version']}")
    print("Files:", clinvar_data["files"])

    # 3. Apply a transformation to gnomAD
    print("\nApplying normalization to gnomAD dataset...")
    transformed = dm.transform_dataset("gnomAD", transformation="normalize")
    print("Transformed files:", transformed["files"])

    print("\nDemo complete.")

if __name__ == "__main__":
    main()

