import asyncio
from unittest.mock import patch, AsyncMock
from omnibioai.services.agentic_ai_service import AgenticAIService
from omnibioai.services.logger_service import logger

# ----------------------------
# Test dataset metadata
# ----------------------------
dataset_metadata = {"id": "ds001", "task_type": "single_cell"}
plugin_outputs = ["Clustering"]

# ----------------------------
# Fully mock RAGService, LLM, ModelZoo, and ExperimentTracking
# ----------------------------
with patch("omnibioai.services.agentic_ai_service.RAGService") as MockRAG, \
     patch("omnibioai.services.agentic_ai_service.ChatOllama") as MockLLM, \
     patch("omnibioai.services.agentic_ai_service.ModelZooService") as MockModelZoo, \
     patch("omnibioai.services.agentic_ai_service.ExperimentTrackingService") as MockTracker, \
     patch("omnibioai.services.agentic_ai_service.initialize_agent") as MockAgentInit:

    # Mock async methods
    mock_rag_instance = MockRAG.return_value
    mock_rag_instance.retrieve_async = AsyncMock(return_value=["Abstract1", "Abstract2"])

    mock_llm_instance = MockLLM.return_value
    mock_llm_instance.generate_async = AsyncMock(return_value="GeneA, GeneB")

    mock_modelzoo_instance = MockModelZoo.return_value
    mock_modelzoo_instance.list_models.return_value = ["Scanpy-ML", "scVI"]

    mock_tracker_instance = MockTracker.return_value
    mock_tracker_instance.log_agentic_ai_suggestions = AsyncMock()

    # Mock the agent
    mock_agent_instance = AsyncMock()
    MockAgentInit.return_value = mock_agent_instance
    mock_agent_instance.arun.return_value = "Suggested next steps"

    # ----------------------------
    # Initialize AgenticAIService
    # ----------------------------
    service = AgenticAIService()

    # ----------------------------
    # Run agentic pipeline
    # ----------------------------
    result = asyncio.run(service.run_agentic_pipeline(dataset_metadata, plugin_outputs, query="BRCA1"))

    # ----------------------------
    # Verify results
    # ----------------------------
    assert "agent_response" in result
    assert "models" in result
    logger.info(f"AgenticAIService returned: {result}")

