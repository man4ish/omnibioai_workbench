# omnibioai/services/agentic_ai_service.py

import asyncio
from functools import lru_cache
from typing import List, Dict
import os
from .logger_service import logger
from .rag_service import RAGService
from .llm_service import LLMService
from .experiment_tracking_service import ExperimentTrackingService
from .model_zoo_service import ModelZooService


class AgenticAIService:
    """
    Agentic AI Service: Provides intelligent recommendations for analysis, 
    experiment planning, gene/pathway exploration, and literature suggestions.
    Integrates with RAG, LLMs, Experiment Tracking, and Model Zoo.
    """

    def __init__(self):
        self.rag = RAGService()
        self.llm = LLMService()
        self.experiment_tracker = ExperimentTrackingService()
        self.model_zoo = ModelZooService()

    async def suggest_next_analysis(self, dataset_metadata: Dict, plugin_outputs: List[str]) -> List[str]:
        """
        Suggest next analysis steps based on dataset type and plugin outputs.
        Returns a list of recommended plugin names.
        """
        suggestions = []
        data_type = dataset_metadata.get("type", "").lower()

        # Rule-based suggestions
        if data_type == "single_cell":
            if "clustering" not in plugin_outputs:
                suggestions.append("Run clustering plugin")
            if "pathway_enrichment" not in plugin_outputs:
                suggestions.append("Run pathway enrichment plugin")
            if "gene_annotation" not in plugin_outputs:
                suggestions.append("Run gene annotation plugin")

        elif data_type == "variant":
            if "variant_annotation" not in plugin_outputs:
                suggestions.append("Run variant annotation plugin")
            if "network_analysis" not in plugin_outputs:
                suggestions.append("Run network analysis plugin")

        elif data_type == "multi_omics":
            suggestions.append("Run ML predictor plugin")
            suggestions.append("Run pathway enrichment plugin")

        logger.info(f"[AgenticAI] Suggested analyses: {suggestions}")
        return suggestions

    async def suggest_genes_or_pathways(self, query: str) -> List[str]:
        """
        Suggest relevant genes or pathways using LLM + RAG.
        """
        abstracts = await self.rag.retrieve_async(query)
        summary = await self.llm.summarize_async(abstracts)
        suggestions = self.extract_genes_pathways(summary)
        logger.info(f"[AgenticAI] Genes/Pathways suggestions for '{query}': {suggestions}")
        return suggestions

    def extract_genes_pathways(self, summary_text: str) -> List[str]:
        """
        Extract genes or pathways from LLM summary using NLP/entity recognition.
        Placeholder: Replace with BioBERT or regex-based extraction.
        """
        # Example static extraction (replace with NLP model in production)
        return ["GeneA", "GeneB", "PathwayX"]

    async def suggest_literature(self, query: str, top_k: int = 5) -> List[Dict]:
        """
        Retrieve top relevant PubMed abstracts or papers.
        """
        abstracts = await self.rag.retrieve_async(query, top_k=top_k)
        logger.info(f"[AgenticAI] Top {top_k} literature suggestions for '{query}' retrieved")
        return abstracts

    @lru_cache(maxsize=128)
    def get_model_recommendations(self, task_type: str) -> List[str]:
        """
        Recommend pre-trained models from Model Zoo for a given task.
        Caching improves response time for repeated queries.
        """
        models = self.model_zoo.list_models(task_type)
        logger.info(f"[AgenticAI] Model recommendations for '{task_type}': {models}")
        return models

    async def suggest_pipeline(self, dataset_metadata: Dict, plugin_outputs: List[str], query: str = "") -> Dict:
        """
        Aggregate suggestions for analyses, genes/pathways, literature, and model recommendations.
        Returns a structured dictionary.
        """
        analyses = await self.suggest_next_analysis(dataset_metadata, plugin_outputs)
        genes_pathways = await self.suggest_genes_or_pathways(query) if query else []
        literature = await self.suggest_literature(query) if query else []
        models = self.get_model_recommendations(dataset_metadata.get("task_type", "general"))

        suggestion_package = {
            "analyses": analyses,
            "genes_pathways": genes_pathways,
            "literature": literature,
            "models": models
        }

        # Track suggestions in Experiment Tracking
        self.experiment_tracker.log_agentic_ai_suggestions(dataset_metadata, suggestion_package)
        return suggestion_package

