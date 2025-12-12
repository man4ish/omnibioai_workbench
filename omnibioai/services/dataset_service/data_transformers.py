from typing import List
import os

class DataTransformers:
    """
    Apply transformations to datasets, e.g., normalization, indexing, annotation prep.
    """

    def apply(self, files: List[str], transformation: str = "normalize") -> List[str]:
        transformed_files = []
        for f in files:
            # Placeholder: real transformations can be added here
            transformed_file = f"{os.path.splitext(f)[0]}_{transformation}{os.path.splitext(f)[1]}"
            # For demo, just copy file
            import shutil
            shutil.copy(f, transformed_file)
            transformed_files.append(transformed_file)
        return transformed_files
