from .client import RAGBioClient
from .loaders import RAGBioLoader
from .exceptions import RAGBioPipelineError

class GeneDiscoveryPipeline:
    """
    Full RAG pipeline:
    - PubMed Search
    - Embedding creation
    - Index building
    - Output loading
    """

    def __init__(self, base_dir):
        self.client = RAGBioClient(base_dir)
        self.loader = RAGBioLoader(base_dir)

    def run(self, query, query_name):
        try:
            # Step 1: Search
            self.client.search_pubmed(query, retmax=500, query_name=query_name)

            # Step 2: Embeddings
            self.client.run_embeddings(query_name)

            # Step 3: Build index
            self.client.build_index(query_name)

            # Step 4: Load data
            abstracts = self.loader.load_abstracts(query_name)
            embeddings = self.loader.load_embeddings(query_name)

            return {
                "abstracts": abstracts,
                "embeddings": embeddings,
            }

        except Exception as e:
            raise RAGBioPipelineError(str(e))
