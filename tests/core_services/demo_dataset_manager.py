"""
Demo for DatasetManager
"""

from omnibioai.services.dataset_service.dataset_manager import DatasetManager
import os

def main():
    # Ensure a demo data folder exists
    demo_data_dir = "data"
    os.makedirs(demo_data_dir, exist_ok=True)
    
    # Create dummy files for the demo
    sample_files = [os.path.join(demo_data_dir, f"sample{i}.csv") for i in range(1, 3)]
    for f in sample_files:
        with open(f, "w") as fh:
            fh.write("col1,col2\n1,2\n3,4\n")
    
    # Initialize DatasetManager
    manager = DatasetManager(cache_dir="data/cache")

    # Register a local dataset (optional)
    manager.register_dataset(name="demo_dataset", files=sample_files, version="v1.0")

    # Retrieve the dataset
    dataset_info = manager.get_dataset("demo_dataset")
    print("Original dataset info:")
    print(dataset_info)

    # Transform the dataset
    transformed_info = manager.transform_dataset("demo_dataset", transformation="normalize")
    print("Transformed dataset info:")
    print(transformed_info)

if __name__ == "__main__":
    main()

