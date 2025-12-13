"""
Module: rag_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides a wrapper around the RAGAssistant pipeline for core retrieval-augmented generation (RAG) services.
    Simplifies querying and retrieving summaries and citations for questions.

Usage:
    from omnibioai.services.rag_service import RAGServiceCore

    rag_service = RAGServiceCore()
    result = rag_service.query("What are the key genes involved in metabolism?")
    print(result["summary"])
    print(result["citations"])

Classes:
    - RAGServiceCore:
        Core service class wrapping RAGAssistant for simplified access.
        
        Methods:
            * __init__():
                Initializes the RAGAssistant instance.
            * query(question: str) -> dict:
                Runs the RAGAssistant on the given question and returns a dictionary
                containing 'summary' and 'citations'.

Dependencies:
    - ragbio.pipeline.rag_pipeline.RAGAssistant: Core RAG pipeline implementation.
"""

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
