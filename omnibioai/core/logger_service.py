"""
logging_setup.py

Configures centralized logging for the OmniBioAI application.

This module sets up a logger with the following features:
- Logs messages to both a file and the console (stream handler).
- Log file location is defined by LOG_DIR from the configuration.
- Logging level is set via LOG_LEVEL from the configuration.
- Standard log format: timestamp, log level, logger name, and message.
- Provides a reusable logger instance named 'OmniBioAI'.

Usage:
    from .logging_setup import logger
    logger.info("This is an info message")
    logger.error("This is an error message")
"""

import logging
from .config import LOG_LEVEL, LOG_DIR
import os

LOG_FILE = os.path.join(LOG_DIR, "omnibioai.log")

logging.basicConfig(
    level=LOG_LEVEL,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger("OmniBioAI")
