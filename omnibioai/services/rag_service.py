from ragbio.pipeline.rag_pipeline import RAGAssistant

class RAGServiceCore:
    """Wrapper around RAGAssistant for core service."""

    def __init__(self):
        self.assistant = RAGAssistant()

    def query(self, question):
        result = self.assistant.run(question)  # or the appropriate method
        return {
            "summary": result.get("summary", ""),
            "citations": result.get("citations", [])
        }
