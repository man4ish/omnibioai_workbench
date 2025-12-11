from .client import RAGBioClient

class RAGBioEmbeddingService:
    """
    Wraps embedding generation and retrieval.
    """

    def __init__(self, base_dir):
        self.client = RAGBioClient(base_dir)

    def generate(self, query_name):
        return self.client.run_embeddings(query_name)
