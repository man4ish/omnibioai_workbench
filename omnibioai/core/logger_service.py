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
