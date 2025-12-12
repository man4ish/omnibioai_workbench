from omnibioai.services.logger_service import logger
from omnibioai.services.llm_service import LLMService

# Initialize service
llm_service = LLMService()

# Test prompt
prompt_text = "Summarize the gene BRCA1 in one sentence."
response = llm_service.prompt(prompt_text)

logger.info(f"LLM response: {response}")
print("LLM response:", response)

