from omnibioai.services.dataset_service.data_transformers import DataTransformers
import os

def main():
    # Create dummy files
    os.makedirs("data", exist_ok=True)
    files = ["data/sample1.csv", "data/sample2.csv"]
    for f in files:
        with open(f, "w") as fh:
            fh.write("col1,col2\n1,2\n3,4\n")

    # Apply transformations
    transformer = DataTransformers()
    transformed_files = transformer.apply(files, transformation="normalize")
    print(transformed_files)

if __name__ == "__main__":
    main()
