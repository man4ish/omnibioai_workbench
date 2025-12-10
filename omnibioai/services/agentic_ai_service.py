# omnibioai/services/agentic_ai_service.py

import asyncio
from functools import lru_cache
from typing import List, Dict
import os
from .logger_service import logger
from .rag_service import RAGService
from .experiment_tracking_service import ExperimentTrackingService
from .model_zoo_service import ModelZooService

# LangChain & Ollama imports
from langchain.chat_models import ChatOllama
from langchain.schema import HumanMessage
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate

# LangGraph imports
from langgraph import Node, Graph


class AgenticAIService:
    """
    Agentic AI Service: Provides intelligent recommendations for analysis, 
    experiment planning, gene/pathway exploration, and literature suggestions.
    Integrates with RAG, LLMs, Experiment Tracking, Model Zoo, LangChain, and LangGraph.
    """

    def __init__(self):
        # Core services
        self.rag = RAGService()
        self.experiment_tracker = ExperimentTrackingService()
        self.model_zoo = ModelZooService()

        # Ollama LLM via LangChain
        self.llm = ChatOllama(model="deepseek-r1")  

        # LangChain prompt template for analysis suggestions
        self.analysis_prompt = PromptTemplate(
            input_variables=["metadata", "plugin_outputs", "candidate_steps"],
            template=(
                "Dataset metadata: {metadata}\n"
                "Current plugin outputs: {plugin_outputs}\n"
                "Candidate next steps: {candidate_steps}\n"
                "Suggest which analyses to run next and why."
            )
        )
        self.analysis_chain = LLMChain(llm=self.llm, prompt=self.analysis_prompt)

        # LangGraph for plugin dependency management
        self.plugin_graph = Graph()
        self._init_plugin_graph()

    def _init_plugin_graph(self):
        """Initialize LangGraph with plugin dependencies."""
        clustering = Node("Clustering")
        pathway = Node("Pathway Enrichment")
        gene_annotation = Node("Gene Annotation")
        variant_annotation = Node("Variant Annotation")
        network_analysis = Node("Network Analysis")
        ml_predictor = Node("ML Predictor")

        # Example dependencies
        self.plugin_graph.add_node(clustering)
        self.plugin_graph.add_node(pathway)
        self.plugin_graph.add_node(gene_annotation)
        self.plugin_graph.add_node(variant_annotation)
        self.plugin_graph.add_node(network_analysis)
        self.plugin_graph.add_node(ml_predictor)

        self.plugin_graph.add_edge(clustering, pathway)  # pathway after clustering
        self.plugin_graph.add_edge(gene_annotation, pathway)  # pathway needs gene annotation

    async def suggest_next_analysis(self, dataset_metadata: Dict, plugin_outputs: List[str]) -> List[str]:
        """
        Suggest next analysis steps using:
        1) LangGraph dependency management (rule-based)
        2) LLM refinement via LangChain + Ollama
        """
        # ----- Step 1: LangGraph candidate steps -----
        candidate_steps = self.plugin_graph.get_ready_nodes(completed=plugin_outputs)
        logger.info(f"[AgenticAI] Candidate steps from graph: {candidate_steps}")

        if not candidate_steps:
            return []

        # ----- Step 2: LLM refinement -----
        llm_input = {
            "metadata": dataset_metadata,
            "plugin_outputs": plugin_outputs,
            "candidate_steps": candidate_steps
        }
        refined_suggestions = await self.analysis_chain.arun(llm_input)
        # Split into list if returned as comma-separated string
        if isinstance(refined_suggestions, str):
            refined_suggestions = [s.strip() for s in refined_suggestions.split(",")]

        # Remove duplicates and log
        suggestions = list(dict.fromkeys(candidate_steps + refined_suggestions))
        logger.info(f"[AgenticAI] Combined suggestions (LangGraph + LLM): {suggestions}")

        return suggestions

    async def suggest_genes_or_pathways(self, query: str) -> List[str]:
        """
        Suggest relevant genes or pathways using RAG + LLM.
        """
        abstracts = await self.rag.retrieve_async(query)
        summary = await self.llm.generate_async(
            f"Summarize these abstracts and extract relevant genes or pathways:\n{abstracts}"
        )
        suggestions = self.extract_genes_pathways(summary)
        logger.info(f"[AgenticAI] Genes/Pathways suggestions for '{query}': {suggestions}")
        return suggestions

    def extract_genes_pathways(self, summary_text: str) -> List[str]:
        """
        Extract genes or pathways from LLM summary.
        Placeholder: replace with BioBERT/NLP in production.
        """
        return ["GeneA", "GeneB", "PathwayX"]

    async def suggest_literature(self, query: str, top_k: int = 5) -> List[Dict]:
        """
        Retrieve top relevant literature using RAG.
        """
        abstracts = await self.rag.retrieve_async(query, top_k=top_k)
        logger.info(f"[AgenticAI] Top {top_k} literature suggestions for '{query}' retrieved")
        return abstracts

    @lru_cache(maxsize=128)
    def get_model_recommendations(self, task_type: str) -> List[str]:
        """
        Recommend pre-trained models from Model Zoo.
        """
        models = self.model_zoo.list_models(task_type)
        logger.info(f"[AgenticAI] Model recommendations for '{task_type}': {models}")
        return models

    async def suggest_pipeline(self, dataset_metadata: Dict, plugin_outputs: List[str], query: str = "") -> Dict:
        """
        Aggregate suggestions for analyses, genes/pathways, literature, and model recommendations.
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

        # Log suggestions
        self.experiment_tracker.log_agentic_ai_suggestions(dataset_metadata, suggestion_package)
        return suggestion_package
