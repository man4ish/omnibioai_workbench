import json
from pathlib import Path

class DatasetManager:
    """
    Manages local reference datasets with versioning.
    """
    def __init__(self, registry_path="data/reference/dataset_registry.json"):
        self.registry_path = Path(registry_path)
        if not self.registry_path.exists():
            raise FileNotFoundError(f"Dataset registry not found at {self.registry_path}")
        with open(self.registry_path) as f:
            self.registry = json.load(f)

    def get_dataset(self, name: str, version: str = None):
        datasets = self.registry.get(name, [])
        if not datasets:
            raise ValueError(f"Dataset '{name}' not found")
        if version:
            for d in datasets:
                if d["version"] == version:
                    return d
            raise ValueError(f"Version '{version}' of dataset '{name}' not found")
        return datasets[-1]  # return latest by default
