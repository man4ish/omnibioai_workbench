"""
omnibioai.services.rag_service

RAGService: A modular service for OmniBioAI integrating the ragbio package
for biomedical literature exploration using Retrieval-Augmented Generation (RAG).

This module provides a high-level interface to:

1. Download PubMed abstracts using ragbio's data loader.
2. Generate embeddings for abstracts with Ollama and build FAISS similarity indices.
3. Run the RAG pipeline to retrieve relevant abstracts, summarize information, 
   and optionally extract structured drug-target-disease relationships.

The RAGService class exposes three main methods:
    - download_data: Fetches PubMed abstracts for a given search query.
    - generate_embeddings: Creates vector embeddings and FAISS indices for retrieval.
    - run_rag_query: Executes the RAG pipeline to obtain summaries, top PMIDs, 
                     and structured results for a query.

All methods are parameterized for flexibility, allowing customization of:
    - Query name (case study)
    - Search query string
    - Abstract and index storage folders
    - Embedding models
    - Number of top abstracts to retrieve
    - Structured extraction options
    - Output directory for results

Typical usage example:

    from omnibioai.services.rag_service import RAGService

    service = RAGService(base_output_dir="omnibioai_output")
    
    # Step 1: Download PubMed abstracts
    service.download_data(
        query="Alzheimer Disease AND therapy",
        query_name="Alzheimer_CaseStudy"
    )

    # Step 2: Generate embeddings and build FAISS index
    service.generate_embeddings(query_name="Alzheimer_CaseStudy")

    # Step 3: Run RAG query with structured extraction
    result = service.run_rag_query(
        query="Alzheimer Disease AND therapy",
        query_name="Alzheimer_CaseStudy",
        top_k=10,
        structured=True
    )

    print(result["summary"])
    print(result["top_pmids"])
    print(result["structured_results"])

Dependencies:
    - ragbio (for data loading, embeddings, and RAG pipeline)
    - Ollama (embedding and LLM models)
    - FAISS (for similarity search)
    - Python standard libraries: os, subprocess, typing
"""

import os
import subprocess
from typing import Dict, Any
from ragbio.pipeline.rag_pipeline import RAGAssistant
from ragbio.embeddings.embedding_engine import generate_embeddings, build_faiss_index, save_index, load_abstracts
from ragbio import config


class RAGService:
    """
    OmniBioAI RAG Service using ragbio.
    Three modular functions:
    1. download_data
    2. generate_embeddings
    3. run_rag_query
    Each function accepts parameters for flexible use.
    """

    def __init__(self, base_output_dir: str = "output"):
        self.base_output_dir = base_output_dir
        os.makedirs(self.base_output_dir, exist_ok=True)

    # -----------------------------
    # 1️⃣ Download PubMed abstracts
    # -----------------------------
    def download_data(
        self,
        query: str,
        query_name: str = "default",
        retmax: int = 500,
        retstart: int = 0,
        abstract_folder: str = None
    ):
        """
        Download PubMed abstracts using ragbio CLI via subprocess.

        Args:
            query (str): PubMed search query.
            query_name (str): Case study name.
            retmax (int): Maximum abstracts to fetch.
            retstart (int): Starting index for fetching abstracts.
            abstract_folder (str, optional): Folder to store abstracts.
        """
        folder = abstract_folder or config.ABSTRACT_FOLDER
        os.makedirs(folder, exist_ok=True)
        print(f"[INFO] Downloading PubMed abstracts for '{query}' into {folder}...")

        # Call rag_data_loader CLI
        cmd = [
            "python", "-m", "ragbio.utils.rag_data_loader",
            "--search", query,
            "--retmax", str(retmax),
            "--retstart", str(retstart),
            "--query_name", query_name
        ]
        subprocess.run(cmd, check=True)
        print(f"[OK] PubMed abstracts downloaded for query_name={query_name}")

    # -----------------------------
    # 2️⃣ Generate embeddings and FAISS index
    # -----------------------------
    def generate_embeddings(
        self,
        query_name: str,
        abstract_folder: str = None,
        index_folder: str = None,
        model_name: str = None
    ):
        """
        Generate embeddings and build FAISS index for a case study.

        Args:
            query_name (str): Case study name.
            abstract_folder (str, optional): Folder containing abstracts.
            index_folder (str, optional): Folder to store FAISS index and PMID map.
            model_name (str, optional): Ollama embedding model.
        """
        abstracts_folder = abstract_folder or config.ABSTRACT_FOLDER
        idx_folder = index_folder or config.INDEX_FOLDER
        model = model_name or config.MODEL_NAME

        print(f"[INFO] Loading abstracts from {abstracts_folder}/{query_name}...")
        abstracts, pmids = load_abstracts(query_name=query_name)
        if not abstracts:
            raise ValueError(f"No abstracts found for query_name={query_name}")

        print(f"[INFO] Generating embeddings using model {model}...")
        checkpoint_path = os.path.join(idx_folder, query_name, "embedding_checkpoint.json")
        embeddings = generate_embeddings(
            texts=abstracts,
            model_name=model,
            checkpoint_path=checkpoint_path
        )

        print(f"[INFO] Building FAISS index...")
        index = build_faiss_index(embeddings)

        print(f"[INFO] Saving FAISS index and PMID map...")
        save_index(index, pmids, query_name=query_name)
        print(f"[OK] Embeddings and index ready for query_name={query_name}")

    # -----------------------------
    # 3️⃣ Run RAG query
    # -----------------------------
    def run_rag_query(
        self,
        query: str,
        query_name: str = "default",
        top_k: int = 10,
        structured: bool = False,
        output_dir: str = None
    ) -> Dict[str, Any]:
        """
        Run the RAG pipeline for a given query.

        Args:
            query (str): Search query string.
            query_name (str): Case study name.
            top_k (int): Number of top abstracts to retrieve.
            structured (bool): Whether to extract structured drug-target-disease info.
            output_dir (str, optional): Directory to save outputs.

        Returns:
            Dict[str, Any]: Summary, top PMIDs, structured results.
        """
        out_dir = output_dir or os.path.join(self.base_output_dir, query_name)
        os.makedirs(out_dir, exist_ok=True)

        assistant = RAGAssistant(output_dir=out_dir, query_name=query_name)
        summary, top_pmids, structured_results = assistant.run_pipeline(
            query=query,
            top_k=top_k,
            structured=structured
        )

        return {
            "query": query,
            "summary": summary,
            "top_pmids": top_pmids,
            "structured_results": structured_results
        }
