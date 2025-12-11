import json
from pathlib import Path

class RAGBioLoader:
    """
    Loads outputs from RAGBIO queries into Python objects
    for downstream Django services.
    """

    def __init__(self, base_dir):
        self.base_dir = Path(base_dir)

    def load_abstracts(self, query_name):
        abstracts_file = self.base_dir / query_name / "abstracts.json"
        if not abstracts_file.exists():
            return []

        with open(abstracts_file, "r") as f:
            return json.load(f)

    def load_embeddings(self, query_name):
        emb_file = self.base_dir / query_name / "embeddings.json"
        if not emb_file.exists():
            return {}

        with open(emb_file, "r") as f:
            return json.load(f)
