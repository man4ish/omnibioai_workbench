from ragbio.pipeline.rag_pipeline import run_pipeline

class RAGServiceCore:
    """Wrapper around existing ragbio pipeline for core service."""
    
    def query(self, question):
        result = run_pipeline(question)
        return {
            "summary": result.get("summary", ""),
            "citations": result.get("citations", [])
        }

