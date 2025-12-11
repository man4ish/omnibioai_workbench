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
from langchain.agents import Tool, initialize_agent, AgentType
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate

# LangGraph imports
from langgraph import Node, Graph


class AgenticAIService:
    """
    Agentic AI Service with fully asynchronous LangChain agent for
    dynamic, step-by-step experiment planning.
    """

    def __init__(self):
        # Core services
        self.rag = RAGService()
        self.experiment_tracker = ExperimentTrackingService()
        self.model_zoo = ModelZooService()

        # Ollama LLM via LangChain
        self.llm = ChatOllama(model="deepseek-r1")  

        # LangGraph for plugin dependency management
        self.plugin_graph = Graph()
        self._init_plugin_graph()

        # Setup tools for LangChain agent
        self.tools = [
            Tool(
                name="Suggest Analysis",
                func=self._agent_suggest_analysis,
                description="Suggests the next analysis steps based on current plugin outputs and dataset type"
            ),
            Tool(
                name="Suggest Genes/Pathways",
                func=self._agent_suggest_genes_pathways,
                description="Suggests relevant genes or pathways for a given query"
            ),
            Tool(
                name="Suggest Literature",
                func=self._agent_suggest_literature,
                description="Retrieves relevant PubMed abstracts or papers"
            )
        ]

        # Initialize LangChain agent (asynchronous)
        self.agent = initialize_agent(
            tools=self.tools,
            llm=self.llm,
            agent=AgentType.CHAT_CONVERSATIONAL_REACT_DESCRIPTION,
            verbose=True
        )

    def _init_plugin_graph(self):
        """Initialize LangGraph with plugin dependencies."""
        clustering = Node("Clustering")
        pathway = Node("Pathway Enrichment")
        gene_annotation = Node("Gene Annotation")
        variant_annotation = Node("Variant Annotation")
        network_analysis = Node("Network Analysis")
        ml_predictor = Node("ML Predictor")

        self.plugin_graph.add_nodes([clustering, pathway, gene_annotation,
                                     variant_annotation, network_analysis, ml_predictor])

        self.plugin_graph.add_edge(clustering, pathway)
        self.plugin_graph.add_edge(gene_annotation, pathway)

    # ----- Tools for LangChain agent -----
    async def _agent_suggest_analysis(self, input_text: str) -> str:
        """
        Tool for LangChain agent to suggest next plugin steps.
        Expects input_text as JSON-like string with dataset metadata & plugin outputs.
        """
        import json
        try:
            info = json.loads(input_text)
            dataset_metadata = info.get("metadata", {})
            plugin_outputs = info.get("plugin_outputs", [])
        except Exception as e:
            return f"Error parsing input: {e}"

        candidate_steps = self.plugin_graph.get_ready_nodes(completed=plugin_outputs)
        return ", ".join(candidate_steps) if candidate_steps else "No further steps available"

    async def _agent_suggest_genes_pathways(self, query: str) -> str:
        genes = await self.suggest_genes_or_pathways(query)
        return ", ".join(genes)

    async def _agent_suggest_literature(self, query: str) -> str:
        abstracts = await self.suggest_literature(query, top_k=5)
        titles = [a.get("title", "No title") for a in abstracts]
        return ", ".join(titles)

    # ----- Existing methods -----
    async def suggest_genes_or_pathways(self, query: str) -> List[str]:
        abstracts = await self.rag.retrieve_async(query)
        summary = await self.llm.generate_async(
            f"Summarize these abstracts and extract relevant genes or pathways:\n{abstracts}"
        )
        return ["GeneA", "GeneB", "PathwayX"]  # placeholder

    async def suggest_literature(self, query: str, top_k: int = 5) -> List[Dict]:
        return await self.rag.retrieve_async(query, top_k=top_k)

    @lru_cache(maxsize=128)
    def get_model_recommendations(self, task_type: str) -> List[str]:
        return self.model_zoo.list_models(task_type)

    # ----- Agent-driven pipeline -----
    async def run_agentic_pipeline(self, dataset_metadata: Dict, plugin_outputs: List[str], query: str = "") -> Dict:
        """
        Run the LangChain agent to dynamically plan next steps.
        """
        input_prompt = {
            "metadata": dataset_metadata,
            "plugin_outputs": plugin_outputs,
            "query": query
        }

        # Run agent asynchronously
        agent_response = await self.agent.arun(input_prompt)

        # Log agent response
        logger.info(f"[AgenticAI Agent] Response: {agent_response}")

        # Return structured result
        suggestion_package = {
            "agent_response": agent_response,
            "models": self.get_model_recommendations(dataset_metadata.get("task_type", "general"))
        }

        self.experiment_tracker.log_agentic_ai_suggestions(dataset_metadata, suggestion_package)
        return suggestion_package
