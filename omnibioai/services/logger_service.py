"""
Module: logger_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides a centralized logger for the OmniBioAI project using Python's built-in logging module.
    Configures logging to output messages with timestamps and severity levels to the console.
    Also includes a small utility function to ensure directories exist.

Usage:
    from omnibioai.services.logger_service import logger, ensure_dir

    # Logging example
    logger.info("This is an info message.")
    logger.error("This is an error message.")

    # Directory utility example
    ensure_dir("data/outputs")

Attributes:
    - logger: Configured logger instance for OmniBioAI.
    
Functions:
    - ensure_dir(path: str):
        Ensures that the directory at `path` exists. Creates it if missing.
"""


import logging

logger = logging.getLogger("OmniBioAI")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Optional utility functions
def ensure_dir(path):
    import os
    if not os.path.exists(path):
        os.makedirs(path)

