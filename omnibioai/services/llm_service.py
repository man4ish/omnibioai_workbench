import ollama
from omnibioai.core.config import RAG_LLM_MODEL
from omnibioai.services.logger_service import logger


class LLMService:
    """
    Service to interact with Ollama LLMs (DeepSeek / LLaMA3).
    Uses Ollama server running on localhost:11434
    """

    def __init__(self, model_name: str = RAG_LLM_MODEL):
        self.model_name = model_name
        logger.info(f"Initialized LLMService with model: {self.model_name}")

    def prompt(
        self,
        text: str,
        temperature: float = 0.7,
        max_tokens: int = 512,
    ) -> str:
        """
        Send a prompt to Ollama and return generated text.
        """
        try:
            response = ollama.generate(
                model=self.model_name,
                prompt=text,
                options={
                    "temperature": temperature,
                    "num_predict": max_tokens,
                },
            )

            output = response.get("response", "").strip()
            logger.info("LLM response generated successfully")
            return output

        except Exception as e:
            logger.error(f"LLM generation failed: {e}")
            return ""
