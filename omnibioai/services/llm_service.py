"""
Module: llm_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    This module provides the LLMService class for interacting with Ollama 
    large language models (LLMs), such as DeepSeek or LLaMA3. It allows 
    sending prompts to the model and retrieving generated responses.

Usage:
    from llm_service import LLMService

    # Initialize the service
    llm = LLMService(model_name="DeepSeek")

    # Generate a response
    response_text = llm.prompt(
        text="Explain CRISPR gene editing.",
        temperature=0.7,
        max_tokens=512
    )

Classes:
    - LLMService: Encapsulates interaction with Ollama LLMs. Provides a
      method `prompt` to send text prompts and receive generated responses.

Dependencies:
    - ollama: Python client for Ollama LLM server.
    - omnibioai.core.config: For default RAG_LLM_MODEL.
    - omnibioai.services.logger_service: For logging.
"""

import ollama
import asyncio
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

    async def generate_async(self, text: str, temperature: float = 0.7, max_tokens: int = 512) -> str:
        """
        Async wrapper for the synchronous `prompt` method.
        """
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(None, self.prompt, text, temperature, max_tokens)        
