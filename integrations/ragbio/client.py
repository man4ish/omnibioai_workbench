import subprocess
from pathlib import Path
from .exceptions import (
    RAGBioConnectionError,
    RAGBioQueryError,
)

class RAGBioClient:
    """
    Adapter to call RAGBIO search, embeddings, indexing and pipelines.

    Example:
        client = RAGBioClient(base_dir="data/PubMed")
        client.search_pubmed("BRCA1 AND cancer")
    """

    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)

    # -------------------------------
    # PubMed Search
    # -------------------------------
    def search_pubmed(self, query: str, retmax=200, query_name="default"):
        cmd = [
            "python",
            "-m",
            "ragbio.utils.rag_data_loader",
            "--search",
            query,
            "--retmax",
            str(retmax),
            "--retstart",
            "0",
            "--query_name",
            query_name,
        ]

        return self._run_cmd(cmd, "PubMed search")

    # -------------------------------
    # Embeddings
    # -------------------------------
    def run_embeddings(self, query_name: str):
        cmd = [
            "python",
            "-m",
            "ragbio.embeddings.embedding_engine",
            "--query_name",
            query_name,
        ]

        return self._run_cmd(cmd, "Embedding generation")

    # -------------------------------
    # Index Build
    # -------------------------------
    def build_index(self, query_name: str):
        index_path = self.base_dir / query_name / "index"
        index_path.mkdir(exist_ok=True, parents=True)

        cmd = [
            "python",
            "-m",
            "ragbio.index.index_engine",
            "--query_name",
            query_name,
        ]

        return self._run_cmd(cmd, "Index building")

    # -------------------------------
    # Utility: run commands
    # -------------------------------
    def _run_cmd(self, cmd, task_name):
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            raise RAGBioQueryError(
                f"{task_name} failed: {e.stderr}"
            )
        except Exception as e:
            raise RAGBioConnectionError(str(e))
