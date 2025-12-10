from ragbio import run_rag_query

class RAGServiceCore:
    """Wrapper around existing ragbio pipeline for core service."""
    
    def query(self, question):
        result = run_rag_query(question)
        return {
            "summary": result.get("summary", ""),
            "citations": result.get("citations", [])
        }

