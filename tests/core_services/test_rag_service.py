from unittest.mock import patch, MagicMock
from omnibioai.services.rag_service import RAGServiceCore
from omnibioai.services.logger_service import logger

test_question = "What is the function of the BRCA1 gene?"

# Fully mock the RAGAssistant constructor and its run method
with patch("omnibioai.services.rag_service.RAGAssistant") as MockAssistant:
    # Mock instance returned whenever RAGAssistant() is called
    mock_instance = MockAssistant.return_value

    # Mock run() to return a fake summary
    mock_instance.run.return_value = {
        "summary": "BRCA1 repairs DNA.",
        "citations": ["PMID12345"]
    }

    # Initialize RAGServiceCore (constructor uses the mocked RAGAssistant)
    rag_service = RAGServiceCore()

    # Query using the wrapper
    response = rag_service.query(test_question)

    # Verify the mocked response
    if "summary" in response and "citations" in response:
        logger.info(f"RAGService (mocked) returned keys correctly: {list(response.keys())}")
        logger.info(f"Summary: {response['summary']}")
        logger.info(f"Citations: {response['citations']}")
    else:
        logger.error(f"RAGService response missing expected keys: {response.keys()}")
